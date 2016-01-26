# Invariant checking and validation

Some classes or API specifications may admit a valdidity or consistency/invariant check. For
uniformity in tests, these checks should conform to the API described below.

## Validity checks

Validity tests confirm that members have valid values and that any class invariants have
been preserved.

Classes which offer a validity test should provide a public method

````cpp
X is_valid() const
````

with a return value that converts implicitly to `bool`, giving false if and
only if the validation check fails. The return value may contain additional
information, depending on the type.

## The `check_valid_api` class

The `check_valid_api` template class definied in `check_valid.h` provides additional functionality.
A class `T` deriving from `check_valid<T>` provides additional member functions:

````cpp
void check_valid() const;

void check_valid(const std::string &message) const;

template <typename Exception,typename... Args>
void check_valid_ex(Args &&... args) const;

void assert_valid() const;
````

The first two will throw an exception of type `rdmini::validation_failure` if
`is_valid()` returns a value that implicity converts to false; if no message is
explicitly given, and the return
type of `is_valid` has a member function with signature `R what() const` with return type
`R` convertible to `std::string`, this message is used to construct the
`validation_failure` exception.

The third form will throw an exception of type specified in the template
parameter, constructed from the supplied arguments.

The `assert_valid()` member function will do nothing if `NDEBUG` is defined,
but otherwise will abort with a message if the validation check fails.

A helper class `valid_info` is provided as a convenient return type for
an `is_valid` member function. `valid_info` converts implicitly to and from
`bool`, but also may be constructed from a `std::string` argument describing
how a validation check fails.

## Other functionality in `check_valid.h`

This header also defines two macros, `SOURCE_LINE` (a string literal) and
`SOURCE_LINE_FUNC` (a temporary value of type `const char *`) which expand to the
current source file and line number, or file, line number and function name
respectively.

In addition, two scoped validation check functions are provided, `check_valid_guard`
and `assert_valid_guard`. These produce a guard object (of implementation-specific type)
which will perform validity checks on the supplied object or object pointer at
construction and also when the guard falls out of scope (i.e. is destroyed.)

`check_valid_guard` will throw a `rdmini::validation_failure` exception on check
failure. 

`assert_valid_guard` will perform no checks if `NDEBUG` is defined, but otherwise
will print a message to `stderr` and abort if the check fails.

Example usage:

````cpp
using namespace rdmini;

struct my_class: check_valid_api<my_class> {
    void foo() {
        // abort if is_valid() returns false at beginning or end
        // of the methos.
        auto _(assert_valid_guard(this));

        unsafe_operation();
        // ...
    }

    void unsafe_operation() { /* ... */ }

    valid_info is_valid() const {
        if (a!=b) return "a not equal to b";
        else return true;
    }

    int a,b; // invariant: a==b
};

void bar(my_class &x) {
    // throw if x is invalid at beginning or end of bar()
    auto _(check_valid_guard(x));

    x.unsafe_operation();
}
````

