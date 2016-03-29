#ifndef IOS_UTIL_H_
#define IOS_UTIL_H_

namespace rdmini {

/** Automatically restore iostream state on scope exit */

struct scoped_ios_format {
    std::ios_base &S;
    std::streamsize precision;
    std::ios::fmtflags flags;

    explicit scoped_ios_format(std::ios_base &S_): S(S_) {
        flags=S.flags();
        precision=S.precision();
    }

    ~scoped_ios_format() {
        S.flags(flags);
        S.precision(precision);
    }
};

}

#endif // ndef IOS_UTIL_H_
