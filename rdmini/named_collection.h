#ifndef NAMED_COLLECTION_H_
#define NAMED_COLLECTION_H_

#include <cstddef>
#include <vector>
#include <unordered_map>
#include <string>
#include <stdexcept>

/* Represent a collection of objects with unique keys
 * given by their 'name' field, of type std::string.
 *
 * Entries can be replaced, but not altered in-place.
 */

template <typename T>
struct named_collection {
    typedef T value_type;
    typedef const T &reference;
    typedef const T &const_reference;
    typedef const T *pointer;
    typedef const T *const_pointer;

    typedef typename std::vector<T>::const_iterator iterator;
    typedef iterator const_iterator;
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;

    typedef std::string key_type;

    iterator begin() const { return items.cbegin(); }
    iterator end() const { return items.cend(); }

    iterator cbegin() const { return items.cbegin(); }
    iterator cend() const { return items.cend(); }

    size_t size() const { return items.size(); }
    bool empty() const { return items.empty(); }
    
    reference operator[](size_type n) const { return items[n]; }
    reference at(size_type n) const { return items.at(n); }
    reference front() const { return items.front(); }
    reference back() const { return items.back(); }

    void clear() { items.clear(); keymap.clear(); }
    
    difference_type index(const key_type &n) const {
        auto j=keymap.find(n);
        return j!=keymap.end()?static_cast<difference_type>(j->second):-1;
    }

    iterator find(const key_type &n) const {
        difference_type i=index(n);
        return i>=0?begin()+i:end();
    }

    reference operator[](const key_type &n) { return *find(n); }
    reference at(const key_type &n) {
        auto i=find(n);
        if (i==end()) throw std::out_of_range("no item with key "+n);
        return *i;
    }

    void insert(const value_type &v) {
        const key_type &k=key(v);
        auto i=index(k);
        if (i>=0) items[i]=v;
        else {
            keymap.insert({{k,items.size()}});
            items.push_back(v);
        }
    }

    key_type unique_key(const key_type &k) {
        key_type uk=k;
        int suffix=0;
        while (index(uk)>=0) uk=k+std::to_string(++suffix);
        return uk;
    }

    static const key_type &key(const value_type &v) {
        return v.name;
    }

    std::vector<T> items;
    std::unordered_map<key_type,size_type> keymap;
};

#endif // ndef NAMED_COLLECTION_H_
