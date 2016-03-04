#ifndef SMALL_MAP_H
#define SMALL_MAP_H

/** Vector-backed small map */

namespace rdmini {

template <typename Key,typename Value,class KeyEqual=std::equal_to<Key>,class Allocator=std::allocator<std::pair<Key,Value>>>
struct small_map {
    typedef Key key_type;
    typedef Value mapped_type;
    typedef std::pair<key_type,mapped_type> value_type;
    typedef KeyEqual key_equal;

    typedef const value_type &const_reference;
    typedef const_reference reference;

    typedef Allocator allocator_type;
private:
    typedef std::vector<value_type,allocator_type> store_type;

public:
    typedef typename store_type::size_type size_type;
    typedef typename store_type::difference_type difference_type;

    typedef typename store_type::const_iterator const_iterator;
    typedef const_iterator iterator;

    small_map(const small_map &) =default;
    small_map(small_map &&) =default;

    explicit small_map(const KeyEqual &eq_=KeyEqual(),
        const Allocator &alloc_=Allocator()): eq(eq_),v(alloc_) {}

    explicit small_map(const Allocator &alloc_): v(alloc_) {}

    small_map(const small_map &other,const Allocator &alloc_)
        : v(other.v,alloc_) {}
    small_map(small_map &&other,const Allocator &alloc_)
        : v(std::move(other.v),alloc_) {}

    template <typename I>
    small_map(I b,I e,const KeyEqual &eq_=KeyEqual(),
        const Allocator &alloc_=Allocator()): eq(eq_),v(alloc_) { insert(b,e); }

    small_map(std::initializer_list<value_type> ilist,
        const KeyEqual &eq_=KeyEqual(),const Allocator &alloc_=Allocator())
        : eq(eq_),v(alloc_) { insert(ilist); }

    small_map &operator=(const small_map &) =default;
    small_map &operator=(small_map &&) =default;

    const_iterator begin() const { return cbegin(); }
    const_iterator cbegin() const { return v.cbegin(); }

    const_iterator end() const { return cend(); }
    const_iterator cend() const { return v.cend(); }
    
    bool empty() const { return v.empty(); }
    size_type size() const { return v.size(); }
    size_type max_size() const { return v.max_size(); }

    void clear() { v.clear(); }

    iterator insert(const value_type &value) {
        auto where=find_in_store(value.first);
        if (where==v.end()) {
            v.push_back(value);
            return std::prev(v.end());
        }
        else {
            *where=value;
            return where;
        }
    }
    
    iterator insert(value_type &&value) {
        auto where=find_in_store(value.first);
        if (where==v.end()) {
            v.push_back(std::move(value));
            return std::prev(v.end());
        }
        else {
            *where=std::move(value);
            return where;
        }
    }

    template <typename I>
    void insert(I b,I e) {
        while (b!=e) insert(*b++);
    }
    
    void insert(std::initializer_list<value_type> ilist) {
        for (const auto &x: ilist) insert(x);
    }

    template <typename... Args>
    std::pair<iterator,bool> emplace(Args &&... args) {
        value_type kv(std::forward<Args>(args)...);
        auto where=find(kv.first);
        if (where!=end()) return std::make_pair(where,false);

        v.push_back(std::move(kv));
        return std::make_pair(std::prev(v.end()),true);
    }

    iterator erase(const_iterator pos) {
/* work asround for defect in libstdc++ for gcc version < 5 */
#if defined(__GNUC__) && __GNUC__ < 5
        return v.erase(v.begin()+std::distance(cbegin(),pos));
#else
        return v.erase(pos);
#endif
    }

    size_type erase(const key_type &key) {
        auto where=find(key);
        if (where==end()) return 0;

        erase(where);
        return 1;
    }

    void swap(small_map &other) {
        std::swap(v,other.v);
    }

    size_type count(const key_type &key) const {
        return find(key)!=end();
    }

    iterator find(const key_type &key) const {
        return find_in_store(key);
    }

    mapped_type &operator[](const Key &key) {
        auto where=find_in_store(key);
        if (where!=v.end()) return where->second;
        v.emplace_back(std::piecewise_construct,std::forward_as_tuple(key),std::tuple<>());
        return std::prev(v.end())->second;
    }

    mapped_type &operator[](Key &&key) {
        auto where=find_in_store(key);
        if (where!=v.end()) return where->second;
        v.emplace_back(std::piecewise_construct,std::forward_as_tuple(std::move(key)),std::tuple<>());
        return std::prev(v.end())->second;
    }

    mapped_type &at(const Key &key) {
        auto where=find_in_store(key);
        if (where!=v.end()) return where->second;
        throw std::out_of_range("missing key");
    }

    const mapped_type &at(const Key &key) const {
        auto where=find_in_store(key);
        if (where!=v.end()) return where->second;
        throw std::out_of_range("missing key");
    }

    KeyEqual key_eq() const { return eq; }
    Allocator get_allocator() const { return v.get_allocator(); }

    friend bool operator==(const small_map &a,const small_map &b) {
        if (a.size()!=b.size()) return false;
        auto bend=b.end();
        for (const auto &e: a) {
            auto bi=b.find(e.first);
            if (bi==bend || e.second!=bi->second) return false;
        }

        return true;
    }

    friend bool operator!=(const small_map &a,const small_map &b) {
        return !(a==b);
    }

private:
    store_type v;
    KeyEqual eq;

    typename store_type::const_iterator find_in_store(const key_type &key) const {
        auto b=v.begin();
        auto e=v.end();
        while (b!=e) if (eq(b->first,key)) break; else ++b;
        return b;
    }

    typename store_type::iterator find_in_store(const key_type &key) {
        auto b=v.begin();
        auto e=v.end();
        while (b!=e) if (eq(b->first,key)) break; else ++b;
        return b;
    }
};

} // namespace rdmini;

#endif // SMALL_MAP_H
