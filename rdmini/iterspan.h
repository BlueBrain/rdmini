#ifndef ITERSPAN_H_
#define ITERSPAN_H_

#include <iterator>
#include <utility>
#include <type_traits>

template <typename T>
class iterspan: private std::pair<T,T> {
    typedef typename std::is_same<std::random_access_iterator_tag,typename std::iterator_traits<T>::iterator_category> is_random_access;

public:
    iterspan(const T &b_,const T &e_): std::pair<T,T>(b_,e_) {}

    T begin() const { return this->first; }
    T end() const { return this->second; }

    typedef T iterator;
    typedef typename std::iterator_traits<T>::value_type value_type;
    typedef typename std::iterator_traits<T>::reference reference;
    typedef typename std::iterator_traits<T>::difference_type difference_type;
    typedef typename std::make_unsigned<difference_type>::type size_type;

    template <typename P=is_random_access>
    typename std::enable_if<P::value,reference>
    operator[](size_type i) const { return begin()[static_cast<difference_type>(i)]; }

    template <typename P=is_random_access>
    typename std::enable_if<P::value,size_type>
    size() const {
        auto d=end()-begin();
        return d<0?0:static_cast<size_type>(d);
    }

};

template <typename T,typename = typename std::iterator_traits<T>::iterator_category>
inline iterspan<T> make_span(T b, T e) { return iterspan<T>(b,e); }

template <typename C,typename = typename std::iterator_traits<typename C::iterator>::iterator_category>
inline iterspan<typename C::iterator> make_span(C &c) { return iterspan<typename C::iterator>(c.begin(),c.end()); }

template <typename C,typename = typename std::iterator_traits<typename C::const_iterator>::iterator_category>
inline iterspan<typename C::const_iterator> make_span(const C &c) { return iterspan<typename C::const_iterator>(c.cbegin(),c.cend()); }

#endif // ndef ITERSPAN_H_

