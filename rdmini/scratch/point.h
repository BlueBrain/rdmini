#ifndef POINT_H_
#define POINT_H_

struct point3d: public std::array<double,3> {
    point3d() = default;

    explicit point3d(double x0) {
        data()[0]=x0;
    }

    point3d(double x0, double x1) {
        data()[0]=x0;
        data()[1]=x1;
    }

    point3d(double x0, double x1, double x2) {
        data()[0]=x0;
        data()[1]=x1;
        data()[2]=x2;
    }

    // arithmetic operations

    point3d &operator+=(const point3d &p) {
        (*this)[0]+=p[0];
        (*this)[1]+=p[1];
        (*this)[2]+=p[2];
        return *this;
    }

    point3d &operator-=(const point3d &p) {
        (*this)[0]-=p[0];
        (*this)[1]-=p[1];
        (*this)[2]-=p[2];
        return *this;
    }

    point3d &operator*=(double x) {
        (*this)[0]*=x;
        (*this)[1]*=x;
        (*this)[2]*=x;
        return *this;
    }

    point3d &operator/=(double x) {
        (*this)[0]/=x;
        (*this)[1]/=x;
        (*this)[2]/=x;
        return *this;
    }

    // equality testing

    bool operator==(const point3d &x) const {
        return (*this)[0]==x[0] && (*this)[1]==x[1] && (*this)[2]==x[2];
    }

    bool operator!=(const point3d &x) const {
        return !(*this==x);
    }

    static point3d zero() { return point3d{0,0,0}; }
};

inline point3d operator+(point3d p, const point3d &q) {
    return p+=q;
}

inline point3d operator-(point3d p, const point3d &q) {
    return p-=q;
}

inline point3d operator*(point3d p, double x) {
    return p*=x;
}

inline point3d operator/(point3d p, double x) {
    return p/=x;
}

inline double dot(const point3d &p, const point3d &q) {
    return p[0]*q[0]+p[1]*q[1]+p[2]*q[2];
}

inline point3d cross(const point3d &p, const point3d &q) {
    return point3d{
        p[1]*q[2] - p[2]*q[1],
        p[2]*q[0] - p[0]*q[2],
        p[0]*q[1] - p[1]*q[0]
    };
}

inline double dist2(point3d a, const point3d &b) {
    a -= b;
    return dot(a,a);
}

inline double distance(point3d a, const point3d &b) {
    return std::sqrt(dist2(a,b));
}

struct bbox3d {
    /** Construct empty bounding box. */
    bbox3d(): empty_flag(true) {}
    
    /** Construct bounding box encompassing a single point. */
    explicit bbox3d(const point3d &x):
        x0(x), x1(x), empty_flag(false) {}
    
    /** Construct bounding box representing the interval [x0,x1]
     * with respect to the partial order on points. */
    bbox3d(const point3d &x0_, const point3d &x1_):
        x0(x0_), x1(x1_), empty_flag(!partial_leq(x0_,x1_)) {}
    
    /** Construct bounding box that encompasses the points descibed
     * by given iterator range. */
    template <typename I>
    bbox3d(I b, I e): empty_flag(true) {
        while (b != e) insert(*b++);
    }

    /** Expand bounding box to minimally enclose point p. */
    void insert(const point3d &p) {
        if (empty()) {
            x0 = x1 = p;
            empty_flag = false;
        }
        else { 
            x0 = meet(x0, p);
            x1 = join(x1, p);
        }
    }

    /* Reset bounding box to empty state. */
    void clear() { empty_flag = true; }
    
    /* Return true if bounding box is empty. */
    bool empty() const { return empty_flag; }

    /* Return true if point p lies within bounding box. */
    bool contains(const point3d &p) const {
        return !empty() && partial_leq(x0,p) && partial_leq(p,x1);
    }

    /** Return reference to minimum point of bounding box.
     *
     * Note that value not well-defined if bounding box is empty.
     */
    const point3d &min() const { return x0; }
    
    /** Return reference to maximum point of bounding box.
     *
     * Note that value not well-defined if bounding box is empty.
     */
    const point3d &max() const { return x1; }
    
private:
    point3d x0, x1;
    bool empty_flag;
    
    static point3d meet(const point3d &a, const point3d &b) {
        return point3d{std::min(a[0],b[0]), std::min(a[1],b[1]), std::min(a[2],b[2])};
    }

    static point3d join(const point3d &a, const point3d &b) {
        return point3d{std::max(a[0],b[0]), std::max(a[1],b[1]), std::max(a[2],b[2])};
    }

    static bool partial_leq(const point3d &a, const point3d &b) {
        return a[0]<=b[0] && a[1]<=b[1] && a[2]<=b[2];
    }
};

#endif // POINT_H_
