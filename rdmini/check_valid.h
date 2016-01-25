#include <cassert>
#include <string>
#include <stdexcept>
#include <utility>

namespace rdmini {

struct validation_failure: std::runtime_error {
    validation_failure(const std::string &what): std::runtime_error(what) {}
    validation_failure(): std::runtime_error("validation failure") {}
};

namespace impl {
#ifdef NDEBUG
    template <typename B>
    struct assert_valid_guard_type {
        assert_valid_guard(const B &me_) {}
    };
#else
    template <typename B>
    struct assert_valid_guard_type {
        assert_valid_guard(const B &me_): me(me_) {
            if (me.is_valid()) return;

            std::fputs("validation failure at guard construction\n",stderr);
            std::abort();
        }

        ~assert_valid_guard() {
            if (me.is_valid()) return;

            std::fputs("validation failure at guard destruction\n",stderr);
            std::abort();
        }

    private:
        const B &me;
    };
#endif

    template <typename B>
    struct check_valid_guard_type {
        check_valid_guard_type(const B &me_): me(me_) {
            auto v=me->is_valid();
            if (!v) throw validation_failure();
        }

        ~check_valid_guard_type() noexcept(false) {
            if (!me->is_valid()) throw validation_failure();
        }

    private:
        const B &me;
    };

    template <typename B,typename Exception>
    struct check_valid_guard_ex_type {
        check_valid_guard_ex_type(const B &me_,Exception &&ex_): me(me_), ex(std::move(ex_)) {
            if (!me.is_valid()) throw ex;
        }

        check_valid_guard_ex_type(const B &me_,const Exception &ex_): me(me_), ex(ex_) {
            if (!me.is_valid()) throw ex;
        }

        ~check_valid_guard_ex_type() noexcept(false) {
            if (!me.is_valid()) throw ex;
        }

    private:
        Exception ex;
        const B &me;
    };

    template <typename X,bool has_message=false>
    struct message_or_default {
        std::string message(const X &,const std::string &default_message) {
            return default_message;
        }
    };

    template <typename X,bool has_message=typename std::is_convertible<decltype(declval(X).message()),std::string>::value>
    struct message_or_default... {
        std::string message(const X &,const std::string &default_message) {
            return default_message;
        }
    };
}

template <typename B>
struct check_valid {
    void check_valid() {
        decltype
    }

    void check_valid(const std::string &message) {
        const B *me=static_cast<const B *>(this);
        if (!me->is_valid()) throw validation_failure(message);
    }

    template <typename Exception=validation_failure,typename... Args>
    void check_valid_ex(Args &&... args) {
        const B *me=static_cast<const B *>(this);
        if (!me->is_valid()) throw Exception(std::forward<Args>(args)...);
    }

    // default implementation of is_valid always returns true; override in derived class T
    bool is_valid() const { return true; }
};

template <typename B>
inline assert_valid_guard_type<B> assert_valid_guard(const B &object) {
    return impl::assert_valid_guard_type<B>(object);
}

template <typename B>
inline assert_valid_guard_type<B> assert_valid_guard(const B *pointer) {
    return impl::assert_valid_guard_type<B>(*pointer);
}

template <typename B>
inline check_valid_guard_type<B> check_valid_guard(const B &object) {
    return impl::check_valid_guard_type<B>(object);
}

template <typename B>
inline check_valid_guard_type<B> check_valid_guard(const B *pointer) {
    return impl::check_valid_guard_type<B>(*pointer);
}

template <typename B,typename Exception>
inline check_valid_guard_type<B> check_valid_guard(const B &object,Exception &&ex) {
    return impl::check_valid_guard_ex_type<B,Exception>(object,std::forward<Exception>(ex));
}

template <typename B,typename Exception>
inline check_valid_guard_type<B> check_valid_guard(const B *pointer,Exception &&ex) {
    return impl::check_valid_guard_ex_type<B,Exception>(*pointer,std::forward<Exception>(ex));
}

} // namespace rdmini

#undef STRINGIFY
#undef STRINGIFY_
#undef SOURCE_LOCATION

#define STRINGIFY_(x) #x
#define STRINGIFY(x) STRINGIFY_(x)
#define SOURCE_LINE __FILE__ ":" STRINGIFY(__LINE__)
#define SOURCE_LINE_FUNC ((std::string(__FILE__ ":" STRINGIFY(__LINE__) ": ")+__func__+"(...)").c_str())
