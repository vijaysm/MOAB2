#ifndef MOAB_MB_ITER_HPP
#define MOAB_MB_ITER_HPP

#define IS_BUILDING_MB
#include "Internals.hpp"
#include "moab/Range.hpp"
#include "moab/Core.hpp"
#include <vector>
#include <algorithm>

struct iBase_EntityArrIterator_Private 
{
  protected:
    iBase_EntityType entType;
    iMesh_EntityTopology entTopo;
    EntityHandle entSet;
    int arrSize;
    
  public:
    iBase_EntityArrIterator_Private( iBase_EntityType type,
                                     iMesh_EntityTopology topology,
                                     EntityHandle set,
                                     int array_size )
      : entType(type), entTopo(topology), entSet(set), arrSize(array_size)
      {}

    virtual ~iBase_EntityArrIterator_Private() {}

    int array_size() const { return arrSize; }

  virtual ErrorCode step(int num_steps, bool &at_end)=0;

    // NOTE: input array must be at least arrLen long
    virtual void get_entities( Core* mb, 
                               EntityHandle* array,
                               int& count_out ) = 0;
    
    virtual ErrorCode reset( Interface* mb ) = 0;
    
    class IsType { 
      private: EntityType type;
      public:  IsType( EntityType t ) : type(t) {}
               bool operator()( EntityHandle h )
                 { return TYPE_FROM_HANDLE(h) == type; }
    };

    void remove_type( std::vector<EntityHandle>& vect, EntityType t ) {
      vect.erase( std::remove_if( vect.begin(), vect.end(), IsType(t) ), vect.end() );
    }
    
    void remove_type( Range& range, EntityType t ) {
      std::pair<Range::iterator,Range::iterator> p = range.equal_range(t);
      range.erase( p.first, p.second );
    }
};


template <class Container> class MBIter : public iBase_EntityArrIterator_Private 
{
  protected:
    Container iterData;
    typename Container::const_iterator iterPos;
    
  public:
    MBIter( iBase_EntityType type,
            iMesh_EntityTopology topology,
            EntityHandle set,
            int array_size )
      : iBase_EntityArrIterator_Private( type, topology, set, array_size ),
        iterPos(iterData.end()) {}
    
    ~MBIter() {}

  typename Container::const_iterator position() const {return iterPos;};

  typename Container::const_iterator end() const {return iterData.end();};

  ErrorCode step(int num_steps, bool &at_end) 
      {
        if (0 > num_steps) return MB_FAILURE;
        
        while (num_steps && iterPos != iterData.end()) {
          num_steps--;
          iterPos++;
        }
        if (iterData.end() == iterPos) at_end = true;
        else at_end = false;

        return MB_SUCCESS;
      }
  
    void get_entities( Core* mb, EntityHandle* array, int& count )
    {
      for (count = 0; count < arrSize && iterPos != iterData.end(); ++iterPos)
        if (mb->is_valid(*iterPos))
          array[count++] = *iterPos;
    }
    
    virtual ErrorCode reset( Interface* mb ) {
      ErrorCode result;
      iterData.clear();
      if (entTopo != iMesh_ALL_TOPOLOGIES) {
        if (entTopo == iMesh_SEPTAHEDRON)
          result = MB_SUCCESS;
        else
          result = mb->get_entities_by_type( entSet, mb_topology_table[entTopo], iterData );
      }
      else if (entType != iBase_ALL_TYPES) {
        result = mb->get_entities_by_dimension( entSet, entType, iterData );
        if (entType == iBase_REGION)
          remove_type( iterData, MBKNIFE );
      }
      else {
        result = mb->get_entities_by_handle( entSet, iterData );
        remove_type( iterData, MBENTITYSET );
        remove_type( iterData, MBKNIFE );
      }
      iterPos = iterData.begin();
      return result;
    }
    
};

typedef MBIter< std::vector<EntityHandle> > MBListIter;
typedef MBIter< moab::Range > MBRangeIter;

#endif
