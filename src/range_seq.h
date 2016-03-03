#include <initializer_list>
#include <utility>
#include <iostream>
#include <vector>
#include <algorithm>

template <typename T>
struct range_seq {
    typedef T value_type;

    std::vector<std::pair<value_type,value_type>> entries;

    template <typename In>
    range_seq(In b,In e) {
        std::vector<value_type> values(b,e);
        std::sort(values.begin(),values.end());
        if (b==e) return;

        std::pair<value_type,value_type> entry{*b,*b};
        while (++b!=e) {
            T n=entry.second;
            ++n;

            T x=*b;
            if (x==n) entry.second=n;
            else {
                entries.push_back(std::move(entry));
                entry=std::pair<T,T>(x,x);
            }
        }

        entries.push_back(std::move(entry));
    }

    template <typename U>
    explicit range_seq(std::initializer_list<U> l): range_seq(l.begin(), l.end()) {}

    friend std::ostream &operator<<(std::ostream &O,const range_seq &rs) {
        bool first=true;
        for (const auto &e: rs.entries) {
            if (!first) O << ',';
            else first=false;

            O << e.first;
            if (e.first!=e.second) O << '-' << e.second;
        }
        return O;
    }
};
