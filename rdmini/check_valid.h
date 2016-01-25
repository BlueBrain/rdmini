#include <cassert>
#include <string>
#include <stdexcept>
#include <utility>
#include <type_traits>

namespace rdmini {

struct validation_failure: std::runtime_error {
    validation_failure(const std::string &what): std::runtime_error(what) {}
    validation_failure(): std::runtime_error("validation failure") {}
};

namespace impl {
    /** Helper routine for has_what<X> type trait.
     *
     * Note: uses SFINAE to test for presence of what() method.
     */

    template <typename X>
    constexpr bool has_what_helper(decltype(X().what()) *p) {
        return std::is_convertible<decltype(*p),std::string>::value;
    }

    template <typename X>
    constexpr bool has_what_helper(...) {
        return false;
    }

    /** Type trait: does X have a method what() returning a value convertible to std::string? */

    template <typename X>
    struct has_what {
        static constexpr bool value=has_what_helper<X>(nullptr);
        using type=std::integral_constant<bool,value>;
    };

    /** Helper class for top-level what() function */

    template <typename X,bool v=has_what<X>::value>
    struct what_helper {
        static std::string what(const X &x) { return x.what(); }
    };

    template <typename X>
    struct what_helper<X,false> {
        static std::string what(const X &x) { return ""; }
    };

    /** Return x.what() as std::string if it exists, or the empty string. */

    template <typename X>
    std::string what(const X &x) {
        return what_helper<X>::what(x);
    }

    /** Abort with what() message. */

    template <typename V>
    inline void abort_what(const V &v,const char *prefix=0) {
        std::fprintf(stderr,"%s%s%s\n",
            prefix?prefix:"",
            (prefix && has_what<V>::value)?": ":"",
            what(v).c_str());
        abort();
    }

    /** Construct validation_failure with what() message, or by default ctor. */

    template <typename V>
    inline validation_failure make_validation_failure(const V &v) {
        if (has_what<V>::value)
            return validation_failure(what(v));
        else
            return validation_failure();
    }

    /** Implementation class for assert_valid_guard function */

#ifdef NDEBUG
    template <typename B>
    struct assert_valid_guard_type {
        assert_valid_guard_type(const B &me_) {}
    };
#else
    template <typename B>
    struct assert_valid_guard_type {
        assert_valid_guard_type(const B &me_): me(me_) {
            if (auto v=me.is_valid()) return;
            else abort_what(v,"validation failure at precondition");
        }

        ~assert_valid_guard_type() {
            if (auto v=me.is_valid()) return;
            else abort_what(v,"validation failure at postcondition");
        }

    private:
        const B &me;
    };
#endif

    /** Implementation class for check_valid_guard function */

    template <typename B>
    struct check_valid_guard_type {
        check_valid_guard_type(const B &me_): me(me_) { run(); }
        ~check_valid_guard_type() noexcept(false) { run(); }

    private:
        const B &me;
        void run() {
            if (auto v=me.is_valid()) return;
            else throw make_validation_failure(v);
        }
    };

} // namespace rdmini::impl


/** CRTP static API for extending functionality of classes with is_valid() method.
 *
 * If a class B provides a method
 *     X is_valid() const
 * where X is implicitly convertible to bool, then deriving B from check_valid_api<B>
 * adds additional methods check_valid() and assert_valid().
 *
 * Refer to validation_api.md document for more details.
 */

template <typename B>
struct check_valid_api {
    /** Check is_valid(), throw validation_failure if false.
     *
     * Use the what() message from the return value of is_valid() if it
     * is available.
     */

    void check_valid() {
        if (auto v=me().is_valid()) return;
        else throw impl::make_validation_failure(v);
    }

    /** Check is_valid(), throw validation_failure if false with supplied message.  */

    void check_valid(const std::string &message) {
        if (!me().is_valid()) throw validation_failure(message);
    }

    /** Check is_valid(), throw user-supplied exception constructed from given arguments
     * if the check fails. */

    template <typename Exception=validation_failure,typename... Args>
    void check_valid_ex(Args &&... args) {
        if (!me().is_valid()) throw Exception(std::forward<Args>(args)...);
    }

    /** Check is_valid(), abort with message if false unless NDEBUG defined. */

    void assert_valid() {
#ifdef NDEBUG
#else
        if (auto v=me().is_valid()) return;
        else impl::abort_what(v,"validation failure");
#endif
    }
private:
    const B &me() const { return static_cast<const B &>(*this); }
};

/** Lightweight utility class for use as a return type for is_valid() methods. */

struct valid_info {
    valid_info() {}
    valid_info(bool v): valid(v) {}
    valid_info(const char *what_): valid(false),what_str(what_) {}

    operator bool() const { return valid; }
    std::string what() const { return what_str?what_str:""; }

private:
    bool valid=false;
    const char *what_str=0;
};

/** Construct validation assertion guard from object reference */

template <typename B>
inline typename std::enable_if<!std::is_pointer<B>::value,impl::assert_valid_guard_type<B>>::type
assert_valid_guard(const B &object) {
    return impl::assert_valid_guard_type<B>(object);
}

/** Construct validation assertion guard from object pointer */

template <typename B>
inline impl::assert_valid_guard_type<B>
assert_valid_guard(const B *pointer) {
    return impl::assert_valid_guard_type<B>(*pointer);
}

/** Construct validation exception guard from object reference */

template <typename B>
inline typename std::enable_if<!std::is_pointer<B>::value,impl::check_valid_guard_type<B>>::type
check_valid_guard(const B &object) {
    return impl::check_valid_guard_type<B>(object);
}

/** Construct validation exception guard from object pointer */

template <typename B>
inline impl::check_valid_guard_type<B>
check_valid_guard(const B *pointer) {
    return impl::check_valid_guard_type<B>(*pointer);
}

} // namespace rdmini

#undef STRINGIFY
#undef STRINGIFY_
#undef SOURCE_LOCATION

#define STRINGIFY_(x) #x
#define STRINGIFY(x) STRINGIFY_(x)
#define SOURCE_LINE __FILE__ ":" STRINGIFY(__LINE__)
#define SOURCE_LINE_FUNC ((std::string(__FILE__ ":" STRINGIFY(__LINE__) ": ")+__func__+"(...)").c_str())
