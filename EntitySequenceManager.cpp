
#ifdef WIN32
#ifdef _DEBUG
// turn off warnings that say they debugging identifier has been truncated
// this warning comes up when using some STL containers
#pragma warning(disable : 4786)
#endif
#endif


#include "PolyEntitySequence.hpp"
#include "EntitySequence.hpp"
#include "EntitySequenceManager.hpp"
#include "MBRange.hpp"
#include "ScdElementSeq.hpp"
#include "ScdVertexSeq.hpp"
#include <assert.h>
#include <iostream>
#include <algorithm>


EntitySequenceManager::EntitySequenceManager()
{
  memset(mLastAccessed, 0, MBMAXTYPE*sizeof(void*));
}

EntitySequenceManager::~EntitySequenceManager()
{
  // delete the entity sequences
  delete_all();
}

void EntitySequenceManager::entity_sequence_created(MBEntitySequence* seq)
{
  mSequenceMap[seq->get_type()].insert(
    std::pair<MBEntityHandle, MBEntitySequence*> (seq->get_start_handle(), seq));
}

void EntitySequenceManager::entity_sequence_deleted(MBEntitySequence* seq)
{
  mSequenceMap[seq->get_type()].erase(seq->get_start_handle());
  mPartlyFullSequenceMap[seq->get_type()].erase(seq->get_start_handle());  
}

void EntitySequenceManager::notify_full(MBEntitySequence* seq)
{
  mPartlyFullSequenceMap[seq->get_type()].erase(seq->get_start_handle());
}

void EntitySequenceManager::notify_not_full(MBEntitySequence* seq)
{
  mPartlyFullSequenceMap[seq->get_type()].insert(
    std::pair<MBEntityHandle, MBEntitySequence*>(seq->get_start_handle(), seq));
}

  //! create a structured sequence of vertices or elements
MBErrorCode EntitySequenceManager::create_scd_sequence(const int imin, const int jmin, const int kmin,
                                                        const int imax, const int jmax, const int kmax,
                                                        const MBEntityType type,
                                                        const int hint_start_id,
                                                        MBEntityHandle &start_handle,
                                                        MBEntitySequence *&seq) 
{
  int this_dim = MBCN::Dimension(type);

    // use > instead of != in the following assert to also catch cases where imin > imax, etc.
  assert((this_dim < 3 || kmax > kmin) &&
         (this_dim < 2 || jmax > jmin) &&
         (this_dim < 1 || imax > imin));

    // compute # entities; not as easy as it would appear...
  int num_ent;
  if (MBVERTEX == type)
    num_ent = (imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1);
  else {
    num_ent = (imax-imin) *
      (this_dim >= 2 ? (jmax-jmin) : 1) *
      (this_dim >= 3 ? (kmax-kmin) : 1);
  }
  
    // get a start handle
  start_handle = get_start_handle(hint_start_id, type, num_ent);
  assert(0 != start_handle);
  if (0 == start_handle) return MB_FAILURE;
  
  if (MBVERTEX == type)
      // create a new structured vertex sequence
    seq = new ScdVertexSeq(this, start_handle, 
                           imin, jmin, kmin, imax, jmax, kmax);
  else
    seq = new ScdElementSeq(this, start_handle, 
                           imin, jmin, kmin, imax, jmax, kmax);
    
  if (NULL == seq) return MB_FAILURE;

  else return MB_SUCCESS;
}

  
/*!
  creates an entity sequence based on number of entities and type.
  uses the hint_start as the start id for the entity handles if possible
  returns the actual start handle and the entity sequence pointer
*/
MBErrorCode EntitySequenceManager::create_entity_sequence( MBEntityType type, int num_ent, int num_nodes,
                                                            int hint_start, MBEntityHandle& start_handle, 
                                                            MBEntitySequence*& seq)
{
  start_handle = get_start_handle(hint_start, type, num_ent);
  
  // actually create the sequence
  return private_create_entity_sequence( start_handle, num_ent, num_nodes, true, seq);
}

MBEntityHandle EntitySequenceManager::get_start_handle(int hint_start, MBEntityType type,
                                                        int num_ent) 
{
  // need to find unused space in the MBEntityHandle ID space
  int dum = 0;
  MBEntityHandle start_hint_handle = CREATE_HANDLE(type, hint_start, dum);
 
  // this is the first of this type we are making 
  if(mSequenceMap[type].empty())
    return start_hint_handle;

  // see if we can use the start hint
  bool can_use_start_hint = false;
  std::map<MBEntityHandle, MBEntitySequence*>::iterator iter =
    mSequenceMap[type].lower_bound(start_hint_handle);
  
  if(iter == mSequenceMap[type].begin())
  {
    if(iter->second->get_start_handle() >= (start_hint_handle + num_ent))
    {
      can_use_start_hint = true;
    }
  }
  else if(iter == mSequenceMap[type].end())
  {
    --iter;
    if(iter->second->get_end_handle() < start_hint_handle)
    {
      can_use_start_hint = true;
    }
  } 
  else 
  {
    std::map<MBEntityHandle, MBEntitySequence*>::iterator jter = iter;
    --jter;
    if((jter->second->get_end_handle() < start_hint_handle) && 
        (iter->second->get_start_handle() >= (start_hint_handle+num_ent)) )
    {
      can_use_start_hint = true;
    }
  } 

  // can we use the start hint?
  if(can_use_start_hint)
    return start_hint_handle;
  else
    return mSequenceMap[type].rbegin()->second->get_end_handle() + 1;
}

MBErrorCode EntitySequenceManager::private_create_entity_sequence(MBEntityHandle start,
                                                           int num_ent, int num_nodes,
                                                           bool full,
                                                           MBEntitySequence *& seq)
{
  MBEntityType type = TYPE_FROM_HANDLE(start);
  if(type == MBVERTEX)
    seq = new VertexEntitySequence(this, start, num_ent, full);
  else if(type == MBPOLYGON || type == MBPOLYHEDRON)
    seq = new PolyEntitySequence(this, start, num_ent, num_nodes, full);
  else
    seq = new ElementEntitySequence(this, start, num_ent, num_nodes, full);
  
  return MB_SUCCESS;
}


MBErrorCode EntitySequenceManager::create_vertex(const double coords[3], MBEntityHandle& handle)
{
    
  // see if there is an existing sequence that can take this new vertex
  if(!mPartlyFullSequenceMap[MBVERTEX].empty())
  {
    MBEntitySequence* seq = mPartlyFullSequenceMap[MBVERTEX].begin()->second;
    handle = seq->get_unused_handle();
    
    static_cast<VertexEntitySequence*>(seq)->
      set_coordinates(handle, coords[0], coords[1], coords[2]);

    return MB_SUCCESS;
  }

  // we need to make a new entity sequence
  if(!mSequenceMap[MBVERTEX].empty())
  {
    handle = mSequenceMap[MBVERTEX].rbegin()->second->get_end_handle() + 1;
  }
  else
  {
    int err=0;
    handle = CREATE_HANDLE(MBVERTEX, MB_START_ID, err);
  }

  VertexEntitySequence* seq = new VertexEntitySequence(this, handle, 4096, false);
  handle = seq->get_unused_handle();
  seq->set_coordinates(handle, coords[0], coords[1], coords[2]);
  
  return MB_SUCCESS;

}


MBErrorCode EntitySequenceManager::create_element(MBEntityType type, 
                                                   const MBEntityHandle *conn, 
                                                   const int num_vertices,
                                                   MBEntityHandle& handle)
{
  // see if there is an existing sequence that can take this new element
  std::map<MBEntityHandle, MBEntitySequence*>::iterator iter;
  for(iter = mPartlyFullSequenceMap[type].begin();
      iter != mPartlyFullSequenceMap[type].end();
      ++iter)
  {
    ElementEntitySequence* seq = dynamic_cast<ElementEntitySequence*>(iter->second);
    if(seq->nodes_per_element() == (unsigned int) num_vertices ||
       seq->nodes_per_element() == 0)
    {
      if (MBPOLYGON == type || MBPOLYHEDRON == type) {
        return dynamic_cast<PolyEntitySequence*>(iter->second)->add_entity(conn, num_vertices, handle);
      }
      else {
        handle = seq->get_unused_handle();
        return seq->set_connectivity(handle, conn, num_vertices);
      }
    }

  }

  // we need to make a new entity sequence
  if(!mSequenceMap[type].empty())
  {
    handle = mSequenceMap[type].rbegin()->second->get_end_handle() + 1;
  }
  else
  {
    int err=0;
    handle = CREATE_HANDLE(type, MB_START_ID, err);
  }

  ElementEntitySequence* seq;
  if (MBPOLYGON == type || MBPOLYHEDRON == type)
    seq = new PolyEntitySequence(this, handle, 0, 0, false);
  else
    seq = new ElementEntitySequence(this, handle, 4096, num_vertices, false);

  handle = seq->get_unused_handle();
  MBErrorCode result = seq->set_connectivity(handle, conn, num_vertices);
  
  return result;
}


MBErrorCode EntitySequenceManager::get_entities(MBEntityType type, MBRange &entities) const
{
  
  //index into the static sequence map to get the sequences according to type only 
  std::map<MBEntityHandle, MBEntitySequence*>::const_iterator beg_seq, end_seq;
  beg_seq = mSequenceMap[type].begin();
  end_seq = mSequenceMap[type].end();

  //for each sequence, get all the entity handles it contains

  for (; beg_seq != end_seq; ++beg_seq)
  {
    const MBEntitySequence* tmp_seq = beg_seq->second;
    tmp_seq->get_entities(entities);
  }

  return MB_SUCCESS;
}

MBErrorCode EntitySequenceManager::get_number_entities(MBEntityType type, int& num_entities) 
{
  num_entities = 0;
  //index into the static sequence map to get the sequences according to type only 
  std::map<MBEntityHandle, MBEntitySequence*>::iterator beg_seq, end_seq;
  beg_seq = mSequenceMap[type].begin();
  end_seq = mSequenceMap[type].end();

  //for each sequence, get all the entity handles it contains

  for (; beg_seq != end_seq; ++beg_seq)
  {
    num_entities += beg_seq->second->number_entities();
  }
  return MB_SUCCESS;
}


void EntitySequenceManager::delete_all()
{
  for(MBEntityType i = MBVERTEX; i<MBMAXTYPE; i++)
  {
    while(!mSequenceMap[i].empty())
      delete mSequenceMap[i].begin()->second;
    mSequenceMap[i].clear();
  } 
}


MBErrorCode EntitySequenceManager::find( MBEntityHandle entity_handle,
                                      MBEntitySequence*& sequence ) const
{
  MBEntityType ent_type = TYPE_FROM_HANDLE(entity_handle);

  // check to see if the sequence is cached
  sequence = mLastAccessed[ent_type];
  if (sequence && sequence->get_start_handle() <= entity_handle &&
     sequence->get_end_handle() >= entity_handle)
  {
    return MB_SUCCESS;
  }
  
  sequence = 0;

  if ( !mSequenceMap[ent_type].empty() )
  {
    if (ID_FROM_HANDLE(entity_handle) == 0) 
      return MB_FAILURE;

    // create an iterator to look for which sequence entity handle will be found in
    // using lower bounds function of map
    std::map<MBEntityHandle, MBEntitySequence*>::const_iterator iter =
      mSequenceMap[ent_type].upper_bound(entity_handle);
    --iter;

    // ensure that the entity is bounded by this sequence.  
    // upper_bound will indicate that this sequence can hold this entity
    // but doesn't guarantee that it does!
    MBEntitySequence* seq = iter->second; 
    if ( (entity_handle >= seq->get_start_handle()) &&
         (entity_handle <= seq->get_end_handle()))
    {
      sequence = seq;
      mLastAccessed[ent_type] = seq;
      return MB_SUCCESS; 
    }
    else
      return MB_ENTITY_NOT_FOUND;
  }

  return MB_FAILURE;
}



MBErrorCode EntitySequenceManager::delete_entity( MBEntityHandle entity )
{
  MBEntitySequence* seq;
  find(entity, seq);
  if(seq != NULL)
  {
    seq->free_handle(entity);
    // leave the sequences around for a while if it is empty
    
    return MB_SUCCESS;
  }
  return MB_FAILURE;
}


#ifdef TEST

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <exodusII.h>

int main()
{
  EntitySequenceManager *seq_manager =  new EntitySequenceManager;

  MBEntityHandle start_ids[4] = {1};
  start_ids[1] = static_cast<MBEntityHandle>((rand()/(RAND_MAX+1.0))*1000 + 100);
  start_ids[2] = static_cast<MBEntityHandle>((rand()/(RAND_MAX+1.0))*1000+start_ids[1]);
  start_ids[3] = static_cast<MBEntityHandle>(start_ids[2]+100);
  std::cout << "start ids " << start_ids[0] << " " << start_ids[1] << " " << start_ids[2] << std::endl;

  MBEntitySequence* e_seqs[3];
  seq_manager->create_entity_sequence(start_ids[0], start_ids[1]-start_ids[0], true, e_seqs[0]);
  seq_manager->create_entity_sequence(start_ids[1], start_ids[2]-start_ids[1], true, e_seqs[1]);
  seq_manager->create_entity_sequence(start_ids[2], start_ids[3]-start_ids[2], true, e_seqs[2]);

  MBTag tags[5] = {0,1,2,3,4};
  int junk_data = 100;
  int tmp;

  for(int i=0; i<1000; i++)
  { 
    junk_data = tmp = rand();
    MBEntityHandle ent_hndl = static_cast<MBEntityHandle>((rand()/(RAND_MAX+1.0))*(start_ids[1]-start_ids[0]));
    e_seqs[0]->set_tag_data(ent_hndl, tags[1],&junk_data, 4, 0);  
    junk_data = rand();
    e_seqs[0]->get_tag_data(ent_hndl, tags[1],&junk_data);  
    assert(tmp == junk_data);
  }  

  // clean up
  delete seq_manager;

  seq_manager =  new EntitySequenceManager;

  char filename[100];
  strcpy(filename, "mdbtest.g");

  float exodus_version;
  int CPU_WORD_SIZE = sizeof(double);  // With ExodusII version 2, all floats
  int IO_WORD_SIZE = sizeof(double);   // should be changed to doubles

  int exodusFile = ex_open ( filename, EX_READ, &CPU_WORD_SIZE,
                         &IO_WORD_SIZE, &exodus_version );

  if(exodusFile == -1)
  {
    std::cout << "couldn't open file " << filename << std::endl;
    return 1;
  }

  int numberDimensions = 0;
  int numberNodes_loading = 0;
  int numberElements_loading = 0;
  int numberElementBlocks_loading = 0; 
  int numberNodeSets_loading = 0;
  int numberSideSets_loading = 0;


  char title[MAX_LINE_LENGTH+1];
  int error = ex_get_init ( exodusFile,
                            title,
                            &numberDimensions,
                            &numberNodes_loading,
                            &numberElements_loading,
                            &numberElementBlocks_loading,
                            &numberNodeSets_loading,
                            &numberSideSets_loading);


  // create a sequence to hold the node coordinates
  MBEntitySequence* nodes;
  seq_manager->create_entity_sequence(1, numberNodes_loading+1, true, nodes);
  double *x_array=0;
  double *y_array=0;
  double *z_array=0;
  
  nodes->get_tag_array_ptr(reinterpret_cast<void**>(&x_array), sizeof(double), 0);
  nodes->get_tag_array_ptr(reinterpret_cast<void**>(&y_array), sizeof(double), 1);
  nodes->get_tag_array_ptr(reinterpret_cast<void**>(&z_array), sizeof(double), 2);

  // read in the coordinates
  error = ex_get_coord( exodusFile, x_array, y_array, z_array );
  std::vector<int> block_ids(numberElementBlocks_loading);
  error = ex_get_elem_blk_ids(exodusFile, &block_ids[0]);

  std::vector<MBEntitySequence*> blocks;
  int elementid = 1;

  for(unsigned int num_blocks = 0; num_blocks < block_ids.size(); num_blocks++)
  {
     int block_handle = block_ids[num_blocks];
     int num_elements, num_nodes_per_element, num_attribs;
     char element_type[MAX_STR_LENGTH+1];

     error = ex_get_elem_block ( exodusFile,
                                      block_handle,
                                      element_type,
                                      &num_elements,
                                      &num_nodes_per_element,
                                      &num_attribs );

    
     MBEntitySequence *block; 
     seq_manager->create_entity_sequence( elementid, num_elements, true, block);
     blocks.push_back( block );
     int* conn = 0;
     blocks[num_blocks]->get_tag_array_ptr(reinterpret_cast<void**>(&conn), sizeof(int)*num_nodes_per_element, 0);
     error = ex_get_elem_conn ( exodusFile,
                                 block_handle,
                                 conn);
     

     elementid += num_elements;
  }

  
  int num_nodes;
  for(num_nodes = 0; num_nodes < numberNodes_loading; num_nodes++)
  {
    double x,y,z;
    nodes->get_tag_data(num_nodes, 0, &x);
    nodes->get_tag_data(num_nodes, 1, &y);
    nodes->get_tag_data(num_nodes, 2, &z);

    //std::cout << "node " << num_nodes << " x=" << x << " y=" << y << " z=" << z << std::endl;

  }

  int entity_handle = 0;
  for(unsigned int j=0; j<blocks.size(); j++)
  {
    MBEntitySequence *seq_pointer = blocks[j];

    int number_elements = seq_pointer->get_num_entities();
  
      int conn_array[8] = {0};

    for(int k=0; k<number_elements; k++)
    {
      seq_pointer->get_tag_data(k+entity_handle,0,conn_array);
      //std::cout << "conn is for element " << k+entity_handle << " is ";
      //for(int kk=0; kk<8; kk++)
      //  std::cout<<" "<<conn_array[kk];
      //std::cout<<std::endl;
      
    }

      entity_handle += number_elements;
  }


  ex_close(exodusFile);

  delete seq_manager;

  std::cout << "success" << std::endl;
  
  return 0;
}

#endif


