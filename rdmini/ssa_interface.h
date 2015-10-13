#ifndef SSA_INTERFACE_H_
#define SSA_INTERFACE_H_

/** Public interface to SSA solvers.
 *
 * An implementation class I should inherit publically from ssa_interface<I>
 * and provide the following types and methods:
 *
 * // Implementation-specific index of distinct populations.
 * typedef pop_index;
 *
 * // Convert from population index to a species id, cell id pair.
 * std::pair<size_t,size_t> pop_to_species(pop_index pop_id)
 *
 * // Convert from species_index and cell_index to a population id.
 * pop_index species_to_pop(size_t species_id,size_t cell_id)
 *
 * // Initialise implementation with rd_model description and initial (sim) time.
 * void initialise_impl(const rd_model &, double t0);
 *
 * // Advance simulation state to time t; actual sim time may be greater.
 * // Return actual sim time.
 * double advance_impl(double t_end);
 *
 * // Advance simulation state by minimal step. Return actual sim time.
 * double advance_impl();
 *
 * // Set a population count.
 * void set_count_impl(pop_index pop_id,size_t count);
 *
 * // Query a population count.
 * void get_count_impl(pop_index pop_id);
 *
 *
 * The ssa_interface<I> base class provides the public interface to the
 * implementation methods, together with some convenience methods.
 *
 * // Set a population count.
 * void set_count(size_t species_id,size_t cell_id,size_t count);
 *
 * // Get a population count.
 * size_t get_count(size_t species_id,size_t cell_id);
 *
 * ...
 */
#if 0
// work in progress
//
template <typename Impl>
struct ssa_interface {
    typedef typename Impl::pop_index pop_index;

    struct count_assign_proxy {
        count_assign_proxy(Impl &impl_,pop_index p_id_): impl(impl),p_id(p_id) {}
        count_assign_proxy &operator=(size_t count) {
            impl.set_count(p_id,count);
            return *this;
        }

    private:
        Impl &impl;
        pop_index p_id;
    };

    // forward to implementation:

    pop_index pop_id(size_t species_id,size_t cell_id) const {
        return impl()->pop_index(species_id,cell_id);
    }

    std::pair<size_t,size_t> from_pop_id(pop_index p_id) const {
        return impl()->from_pop_id(p_id);
    }

    void initialise(const rd_model &M,double t0) {
        impl()->initialise(M,t0);
    }
    
    double advance(double t_end) {
        return impl()->advance(t_end);
    }

    double sim_time() const {
        return impl()->sim_time();
    }
    
    // pop count access
    size_t operator[](pop_index p_id) const {
        return impl->get_count(p_id);
    }

    count_assign_proxy operator[](pop_index p_id) {
        return count_assign_proxy(impl(),p_id);
    }

    // common wrappers:
    size_t species_id(pop_index p_id) const {
        return from_pop_id(p_id).first;
    }

    size_t cell_id(pop_index p_id) const {
        return from_pop_id(p_id).second;
    }

    size_t count(size_t species_id,size_t cell_id) const {
        return (*this)[pop_id(species_id,cell_id)];
    }

    count_assign_proxy count(size_t species_id,size_t cell_id) {
        return (*this)[pop_id(species_id,cell_id)];
    }

private:
    size_t n_species=0;
    size_t n_cells=0;
    Impl *impl() { return static_cast<Impl *>(this); }
};

#endif // WIP

#endif SSA_INTERFACE_H_
