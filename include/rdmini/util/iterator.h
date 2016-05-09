#ifndef RDMINI_ITERATOR_H
#define RDMINI_ITERATOR_H

/** Iterator adaptors and utilities */


#include <utility>

namespace rdmini {

/** Random-access input iterator that just counts up.
 *
 * Acts as an adaptor between things that want integer
 * ranges, and things that want iterator ranges.
 */

template <typename int_type>
struct counting_iterator {
    using difference_type=typename std::make_signed<int_type>::type;
    using value_type=int_type;
    using pointer=const value_type *;
    using reference=value_type;
    using iterator_category=std::random_access_iterator_tag;

    value_type i;

    explicit counting_iterator(int_type i_=0): i(i_) {}

    bool operator==(counting_iterator x) const { return i==x.i; }
    bool operator!=(counting_iterator x) const { return i!=x.i; }

    bool operator<=(counting_iterator x) const { return i<=x.i; }
    bool operator>=(counting_iterator x) const { return i>=x.i; }

    bool operator<(counting_iterator x) const { return i<x.i; }
    bool operator>(counting_iterator x) const { return i>x.i; }

    difference_type operator-(counting_iterator x) const {
        return (difference_type)i-(difference_type)x.i;
    }

    counting_iterator &operator+=(difference_type n) { i+=n; return *this; }
    counting_iterator operator+(difference_type n) const { return counting_iterator(i+n); }
    friend counting_iterator operator+(difference_type n,counting_iterator x) {
        return counting_iterator(x.i+n);
    }

    counting_iterator &operator-=(difference_type n) { i-=n; return *this; }
    counting_iterator operator-(difference_type n) const { return counting_iterator(i-n); }
    friend counting_iterator operator-(difference_type n,counting_iterator x) {
        return counting_iterator(x.i-n);
    }

    reference operator[](difference_type n) const { return i+n; }

    reference operator*() const { return i; }

    counting_iterator &operator++() { return ++i,*this; }
    counting_iterator operator++(int) { counting_iterator x(i); return ++i,x; }

    counting_iterator &operator--() { return --i,*this; }
    counting_iterator operator--(int) { counting_iterator x(i); return --i,x; }
};


/** Present a functor/lambda/function object as an output iterator */

template <typename F>
struct functor_iterator_adaptor {
    explicit functor_iterator_adaptor(const F &f_): p(f_) {}
    
    struct proxy {
        F f;

        proxy(const F &f_): f(f_) {}

        template <typename V>
        void operator=(V &&v) { f(std::forward<V>(v)); }
    } p;

    // output iterator interface
    proxy &operator*() { return p; }
    functor_iterator_adaptor &operator++() { return *this; }
    functor_iterator_adaptor &operator++(int) { return *this; }

    typedef void value_type;
    typedef void difference_type;
    typedef void reference;
    typedef void pointer;
    typedef std::output_iterator_tag iterator_category;
};

template <typename F>
functor_iterator_adaptor<F> functor_iterator(const F &f) { return functor_iterator_adaptor<F>(f); }

} // namespace rdmini

#endif // ndef RDMINI_ITERATOR_H
