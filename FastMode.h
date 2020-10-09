#ifndef _H_FAST_MODE_
#define _H_FAST_MODE_

#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>
#include <PhysTools/histogram.h>
#include <PhysTools/likelihood/likelihood.h>

namespace golemfit {
namespace fastmode {
namespace detail {

    
//some useful type traits
template<typename T>
struct remove_reference_wrapper{ using type=T; };

template<typename T>
struct remove_reference_wrapper<std::reference_wrapper<T>>{ using type=T; };

template<typename T>
struct ensure_reference_wrapper{ using type=typename std::reference_wrapper<T>; };                                                                                                                                                

template<typename T>
struct ensure_reference_wrapper<std::reference_wrapper<T>>{
    using type=typename std::reference_wrapper<T>;
};

template<typename T>
struct dereference {
    T& operator()(T& t) {
        return t;
    }
};

template<typename T>
struct dereference<std::reference_wrapper<T>> {
    T & operator()(std::reference_wrapper<T> t) {
        return t.get();
    }
};

template<std::size_t I = 0, typename f, typename... Tp> 
inline typename std::enable_if<I == sizeof...(Tp), void>::type
apply_to_tuple(f& do_thing, std::tuple<Tp...>& t)
{ }

template<std::size_t I = 0, typename f, typename... Tp> 
inline typename std::enable_if<I < sizeof...(Tp), void>::type
apply_to_tuple(f& do_thing, std::tuple<Tp...>& t)                                                                           
{
    do_thing(std::get<I>(t));
    apply_to_tuple<I + 1, f, Tp...>(do_thing, t); 
}

template<typename X>
class P;

template<typename Action>
struct apply_to_bins {
    Action& action;
    apply_to_bins(Action& action):action(action){}
    template<typename HType>
    void operator()(HType& hist) {
        for(auto it=hist.begin(), end=hist.end(); it!=end; ++it) {
            //const auto& itc = static_cast<const phys_tools::likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(*it);
            //action(itc.entries());
            //P<decltype(*it)> itType;
            action(*it);
        }
    }
};

template<typename Action>
struct apply_to_histogram {
    Action& action;
    apply_to_histogram(Action& action):action(action){}
    template<typename T>
    struct apply_general {
        Action& action;
        apply_general(Action& action):action(action){}
        void operator()(T& t) {
            action(t);
        }
    };

    template<typename... Tp>
    struct apply_general<std::tuple<Tp...>> {
        Action& action;
        apply_general(Action& action):action(action){}
        void operator()(std::tuple<Tp...>& t) {
            apply_to_tuple(action, t);
        }
    };

    template<typename T>
    void operator()(T& t) {
        apply_general<T> ag(action);
        ag(t);
    }
    
};

template<typename Action>
struct apply_to_histogram_bins {
    Action& action;
    apply_to_histogram_bins(Action& action):action(action){}

    template<typename T>
    void operator()(T& t) {
        apply_to_bins<Action> atb(action);
        apply_to_histogram<apply_to_bins<Action>> ath(atb);
        ath(t);
    }
};

template<typename binType>
struct clear_bin {
    template<typename T>
    void operator()(T& t) {
        t = binType(0.);
    }
};

template<typename T> struct extract_data_type { using type=typename T::dataType;};
template<typename... Tp> struct extract_data_type<std::tuple<Tp...>> { using type=typename std::tuple_element<0, std::tuple<Tp...>>::type::dataType;};

template<typename T>
struct clear_histogram {
    void operator()(T& t) {
        clear_bin<typename extract_data_type<T>::type> cb;
        apply_to_histogram_bins<decltype(cb)> clear(cb);
        clear(t);
    }
};

template<typename MetaHistType, typename MetaBinner, typename Combiner, typename Event>
struct meta_event_accumulator {
    MetaHistType& meta_hist;
    MetaBinner binner;
    Combiner combiner;

    std::deque<Event> meta_events;
    
    meta_event_accumulator(MetaHistType& mh, MetaBinner b, Combiner c):meta_hist(mh), binner(b), combiner(c){}
    
    template<typename dataType>
    void pops(std::vector<dataType>& events) {
        //Combiner& c = combiner;
        //std::deque<Event>& me = meta_events;
        
        auto combine_and_store = [&](const phys_tools::likelihood::entryStoringBin<std::reference_wrapper<const Event> > m_events){
            if(m_events.entries().size() > 0) {
                Event e = this->combiner(m_events.entries());
                this->meta_events.push_back(e);
            }
        };
        for(const Event& e : events) {
            binner(meta_hist, e);
        }
        apply_to_histogram_bins<decltype(combine_and_store)> athb(combine_and_store);
        athb(meta_hist);
        clear_histogram<MetaHistType> clear;
        clear(meta_hist);
    }

    template<typename T>
    struct ops {
        meta_event_accumulator<MetaHistType, MetaBinner, Combiner,Event>* this_;
        ops(meta_event_accumulator<MetaHistType, MetaBinner, Combiner,Event>* t):this_(t){}
        void operator()(T& t) {
            auto& itc = static_cast<const phys_tools::likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(t);
            auto events = itc.entries();
            //typedef decltype(events) eventsType;
            this_->pops(events);
        }
    };
    
    template<typename dataType>
    struct ops<std::vector<dataType>> {
        meta_event_accumulator<MetaHistType, MetaBinner, Combiner,Event>* this_;
        ops(meta_event_accumulator<MetaHistType, MetaBinner, Combiner,Event>* t):this_(t){}
        void operator()(std::vector<dataType>& events) {
            this_->pops(events);
        }
    };


    template<typename T>
    void operator()(T& t) {
        ops<T>(this)(t);
    }
};


} // namespace detail

template<typename Event, typename MetaHistType, typename MetaBinner, typename Combiner, typename Histogram>
std::deque<Event> get_fastmode_events(MetaHistType& meta_hist, MetaBinner& meta_binner, Combiner&& combiner, Histogram& hist){
    detail::meta_event_accumulator<MetaHistType, MetaBinner, Combiner, Event> fm(meta_hist, meta_binner, combiner);
    detail::apply_to_histogram_bins<decltype(fm)> ath(fm);
    ath(hist);
    return fm.meta_events;
}

} // namespace fastmode
} // namespace golemfit

#endif
