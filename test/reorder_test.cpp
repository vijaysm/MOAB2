#include "moab/Core.hpp"
#include "moab/ReorderTool.hpp"
#include "TestUtil.hpp"

using namespace moab;

// some tag names
const char GLOBAL_ID_NAME[] = "GLOBAL_ID"; /* global ID assigned to each vtx and quad */
const char SET_IDS_NAME[] = "SET_IDS";     /* global IDs of entities in each set */
const char SET_HANDLES_NAME[] = "SET_HANDLES"; /* handles of entities in each set */
const char CONN_IDS_NAME[] = "CONN_IDS"; /* global IDs of vertices in each quad */
const char CONN_NAME[] = "CONN_HANDLES"; /* handles of vertices in each quad */
const char VAR_INTS_NAME[] = "VAR_LEN_INTS"; /* variable length tag on nodes */
const char BIT_NAME[] = "TEST_BIT_TAG";
const int ENTS_PER_SET = 6;
const int BITS_PER_TAG = 2;

Core* mbcore = 0;
Interface* mb = 0;
Tag order_tag = 0;

const size_t INTERVALS = 6;

/* values for variable-length tag data */
void tag_vals_from_gid( int global_id, std::vector<int>& values )
{
  int i = global_id / (INTERVALS+1);
  int j = global_id % (INTERVALS+1);
  int n = global_id % 5 + 1;
  int vals[]= { i, j, n, i+j, j-2*i };
  values.resize(n);
  std::copy( vals, vals+n, values.begin() );
}

unsigned char bits_from_gid( int global_id )
{
  return global_id % (1<<BITS_PER_TAG);
}

unsigned char order_from_gid( int global_id )
{
  return global_id % 3;
}

void coords_from_gid( int global_id, double coords[3] )
{
  int i = global_id / (INTERVALS+1);
  int j = global_id % (INTERVALS+1);
  coords[0] = i;
  coords[1] = j;
  coords[2] = 0.1*(i+j);
}

void build_mesh();
void check_order_by_sets_and_adj();
void call_reorder();
void check_order();
void check_node_coords();
void check_quad_conn();
void check_set_meshset();
void check_list_meshset();
void check_big_meshset();
void check_handle_tag();
void check_varlen_tag();
void check_bit_tag();

int main()
{
    // Define global MOAB instance for use by all tests
  Core mcore;
  mbcore = &mcore;
  mb = &mcore;

  // if this fails, don't bother with anything else
  if (RUN_TEST(build_mesh))
    return 1;
  
    // this test needs be be run before reordering the mesh
  int errors = 0;
  errors += RUN_TEST(check_order_by_sets_and_adj);
  
    // if reorder returned failure, don't bother doing anything else
  int tmp = RUN_TEST(call_reorder);
  if (tmp)
    return tmp+errors;
    
    // test the core stuff
  errors += RUN_TEST(check_order);
  errors += RUN_TEST(check_node_coords);
  errors += RUN_TEST(check_quad_conn);
  errors += RUN_TEST(check_set_meshset);
  errors += RUN_TEST(check_list_meshset);
  errors += RUN_TEST(check_big_meshset);
  errors += RUN_TEST(check_handle_tag);
  errors += RUN_TEST(check_varlen_tag);
  errors += RUN_TEST(check_bit_tag);
  return errors;
}

void build_mesh( )
{
  const unsigned dense = MB_TAG_CREAT|MB_TAG_DENSE;
  const unsigned sparse = MB_TAG_CREAT|MB_TAG_SPARSE;

  ErrorCode rval;
  
    // get/create various tags
  Tag gid;
  rval = mb->tag_get_handle( GLOBAL_ID_NAME, 1, MB_TYPE_INTEGER, gid, dense );
  CHECK_ERR(rval);

  Tag conn_ids;
  rval = mb->tag_get_handle( CONN_IDS_NAME, 4, MB_TYPE_INTEGER, conn_ids, dense );
  CHECK_ERR(rval);
  
  Tag conn_handles;
  rval = mb->tag_get_handle( CONN_NAME, 4, MB_TYPE_HANDLE, conn_handles, dense );
  CHECK_ERR(rval);

  Tag set_ids;
  rval = mb->tag_get_handle( SET_IDS_NAME, ENTS_PER_SET, MB_TYPE_INTEGER, set_ids, sparse );
  CHECK_ERR(rval);

  Tag set_handles;
  rval = mb->tag_get_handle( SET_HANDLES_NAME, ENTS_PER_SET, MB_TYPE_HANDLE, set_handles, sparse );
  CHECK_ERR(rval);
    
  Tag var_data;
  rval = mb->tag_get_handle( VAR_INTS_NAME, 0, MB_TYPE_INTEGER, var_data, dense|MB_TAG_VARLEN );
  CHECK_ERR(rval);
  
  Tag bit_data;
  rval = mb->tag_get_handle( BIT_NAME, BITS_PER_TAG, MB_TYPE_BIT, bit_data, MB_TAG_CREAT );
  CHECK_ERR(rval);
  
  rval = mb->tag_get_handle( "ORDER", 1, MB_TYPE_INTEGER, order_tag, dense );
  CHECK_ERR(rval);
  
    // create and tag vertices
  std::vector<int> values;
  EntityHandle nodes[(INTERVALS+1)*(INTERVALS+1)];
  for (size_t i = 0; i <= INTERVALS; ++i) {
    for (size_t j = 0; j <= INTERVALS; ++j) {
      size_t idx = i*(INTERVALS+1) + j;
      double coords[3];
      coords_from_gid(idx, coords);
      rval = mb->create_vertex( coords, nodes[idx] );
      CHECK_ERR(rval);
      
      int tagval = idx;
      rval = mb->tag_set_data( gid, nodes+idx, 1, &tagval );
      CHECK_ERR(rval);
      
      tag_vals_from_gid( idx, values );
      const void* ptr = &values[0];
      const int size = values.size();
      rval = mb->tag_set_by_ptr( var_data, nodes+idx, 1, &ptr, &size );
      CHECK_ERR(rval);
      
      unsigned char bits = bits_from_gid( idx );
      rval = mb->tag_set_data( bit_data, nodes+idx, 1, &bits );
      CHECK_ERR(rval);
      
      int group = order_from_gid( idx );
      rval = mb->tag_set_data( order_tag, nodes+idx, 1, &group );
      CHECK_ERR(rval);
    }
  }
  
    // create and tag elements
  EntityHandle quads[INTERVALS*INTERVALS];
  for (size_t i = 0; i < INTERVALS; ++i) {
    for (size_t j = 0; j < INTERVALS; ++j) {
      size_t idx = i * INTERVALS + j;
      size_t n0 =  i    * (INTERVALS+1) + j;
      size_t n1 = (i+1) * (INTERVALS+1) + j;
      size_t n2 = (i+1) * (INTERVALS+1) + j + 1;
      size_t n3 =  i    * (INTERVALS+1) + j + 1;
      EntityHandle conn[4] = { nodes[n0], nodes[n1], nodes[n2], nodes[n3] };
      EntityHandle h;
      rval = mb->create_element( MBQUAD, conn, 4, h );
      CHECK_ERR(rval);
      
      int tagval = idx;
      rval = mb->tag_set_data( gid, &h, 1, &tagval );
      CHECK_ERR(rval);
      
      int ids[4] = { n0, n1, n2, n3 };
      rval = mb->tag_set_data( conn_ids, &h, 1, ids );
      CHECK_ERR(rval);
      
      rval = mb->tag_set_data( conn_handles, &h, 1, conn );
      CHECK_ERR(rval);
      
      int group = order_from_gid( idx );
      rval = mb->tag_set_data( order_tag, &h, 1, &group );
      CHECK_ERR(rval);
      
      quads[idx] = h;
    }
  }
  
    // create a few sets
  for (int i = 0; i < 2; ++i) {
    EntityHandle* from = 0;
    size_t count;
    unsigned flag;
    if (i) {
      from = nodes;
      count = (INTERVALS+1)*(INTERVALS+1);
      flag = MESHSET_SET;
    }
    else {
      from = quads;
      count = INTERVALS*INTERVALS;
      flag = MESHSET_ORDERED;
    }
 
    EntityHandle h;
    rval = mb->create_meshset( flag|MESHSET_TRACK_OWNER, h );
    CHECK_ERR(rval);
    
    EntityHandle ents[ENTS_PER_SET];
    int gids[ENTS_PER_SET];
    for (int j = 0; j < ENTS_PER_SET; ++j) {
      int idx = j+2;
      idx = (idx*idx)%count;
      ents[j] = from[idx];
      gids[j] = idx;
    }
     
    rval = mb->add_entities( h, ents, ENTS_PER_SET );
    CHECK_ERR(rval);
    
    rval = mb->tag_set_data( set_ids, &h, 1, gids );
    CHECK_ERR(rval);
    
    rval = mb->tag_set_data( set_handles, &h, 1, ents );
    CHECK_ERR(rval);
  }
    
    // create a set containing all vertices
  EntityHandle allverts;
  rval = mb->create_meshset( MESHSET_SET, allverts );
  CHECK_ERR(rval);
  rval = mb->add_entities( allverts, nodes, (INTERVALS+1)*(INTERVALS+1) );
  CHECK_ERR(rval);
}

void call_reorder()
{
    // do reorder
  ReorderTool tool(mbcore);
  Tag mapping;
  ErrorCode rval = tool.handle_order_from_int_tag( order_tag, -1, mapping );
  CHECK_ERR(rval);
  rval = tool.reorder_entities( mapping );
  CHECK_ERR(rval);
}

void check_order( EntityType type )
{
  ErrorCode rval;
  
  Tag gid;
  rval = mb->tag_get_handle( GLOBAL_ID_NAME, 1, MB_TYPE_INTEGER, gid );
  CHECK_ERR(rval);
  
  Range ents;
  rval = mb->get_entities_by_type( 0, type, ents );
  CHECK_ERR(rval);
  
  std::vector<int> ids(ents.size());
  rval = mb->tag_get_data( gid, ents, &ids[0] );
  CHECK_ERR(rval);
  
  for (size_t i = 1; i < ids.size(); ++i) {
    CHECK( order_from_gid(ids[i-1]) <= order_from_gid( ids[i] ) );
  }
}

void check_order_by_sets_and_adj()
{
  ErrorCode rval;
  
  std::vector<EntityHandle> quads;
  rval = mb->get_entities_by_dimension( 0, 2, quads );
  CHECK_ERR(rval);
  CHECK(!quads.empty());
  
  // group quads by the ordering assigned in build_mesh()
  std::map<int,Range> groups;
  std::vector<int> group_ids(quads.size());
  rval = mb->tag_get_data( order_tag, &quads[0], quads.size(), &group_ids[0] );
  CHECK_ERR(rval);
  for (size_t i = 0; i < quads.size(); ++i)
    groups[group_ids[i]].insert(quads[i]);
  
  // create sets from groups
  Range sets;
  for (std::map<int,Range>::iterator i = groups.begin(); i != groups.end(); ++i) {
    EntityHandle h;
    rval = mb->create_meshset( MESHSET_SET, h );
    CHECK_ERR(rval);
    rval = mb->add_entities( h, i->second );
    CHECK_ERR(rval);
    sets.insert( h );
  }
  
  // Get ordering assigned by set containment
  Tag neworder;
  ReorderTool tool(mbcore);
  rval = tool.handle_order_from_sets_and_adj( sets, neworder );
  CHECK_ERR(rval);
  
  // check that new quad handles are clustered as expected
  std::vector< std::pair<EntityHandle,EntityHandle> > ranges;
  for (std::map<int,Range>::iterator i = groups.begin(); i != groups.end(); ++i) {
    std::vector<EntityHandle> newh(i->second.size());
    rval = mb->tag_get_data( neworder, i->second, &newh[0] );
    CHECK_ERR(rval);
    std::sort(newh.begin(), newh.end());
    CHECK(newh[0] > 0); // zero implies some quad got left out of the reordering
    std::pair<EntityHandle,EntityHandle> p(newh[0], newh[newh.size()-1]);
    ranges.push_back(p);
  }
  std::sort( ranges.begin(), ranges.end() );
  for (size_t i = 1; i < ranges.size(); ++i) {
    CHECK( ranges[i-1].second < ranges[i].first );
  }
  
  // group vertices as we expect handles to be grouped
  std::map< std::vector<int>, Range > vtxgroups;
  Range verts;
  rval = mb->get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
    Range adj;
    rval = mb->get_adjacencies( &*i, 1, 2, false, adj );
    CHECK_ERR(rval);
    std::vector<int> ids(adj.size());
    rval = mb->tag_get_data( order_tag, adj, &ids[0] );
    CHECK_ERR(rval);
    std::sort( ids.begin(), ids.end() );
    ids.erase( std::unique( ids.begin(), ids.end() ), ids.end() );
    vtxgroups[ids].insert( *i );
  }
  
  // check that new vertex handles are clustered as expected
  ranges.clear();
  std::map< std::vector<int>, Range >::iterator j;
  for (j = vtxgroups.begin(); j != vtxgroups.end(); ++j) {
    std::vector<EntityHandle> newh(j->second.size());
    rval = mb->tag_get_data( neworder, j->second, &newh[0] );
    CHECK_ERR(rval);
    std::sort(newh.begin(), newh.end());
    CHECK(newh[0] > 0); // zero implies some quad got left out of the reordering
    std::pair<EntityHandle,EntityHandle> p(newh[0], newh[newh.size()-1]);
    ranges.push_back(p);
  }
  std::sort( ranges.begin(), ranges.end() );
  for (size_t i = 1; i < ranges.size(); ++i) {
    CHECK( ranges[i-1].second < ranges[i].first );
  }
}  
  

void check_order()
{
  check_order( MBVERTEX );
  check_order( MBQUAD );
}


void check_node_coords()
{
  ErrorCode rval;
  
  Tag gid;
  rval = mb->tag_get_handle( GLOBAL_ID_NAME, 1, MB_TYPE_INTEGER, gid );
  CHECK_ERR(rval);
  
  Range verts;
  rval = mb->get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  
  std::vector<int> ids(verts.size());
  rval = mb->tag_get_data( gid, verts, &ids[0] );
  CHECK_ERR(rval);

  std::vector<double> coords(3*verts.size());
  rval = mb->get_coords( verts, &coords[0] );
  CHECK_ERR(rval);
  
  std::vector<double> expected(3*verts.size());
  for (size_t i = 0; i < ids.size(); ++i) 
    coords_from_gid( ids[i], &expected[3*i] );
  
  CHECK_EQUAL( expected, coords );
}


void check_quad_conn()
{
  ErrorCode rval;
  
  Tag gid;
  rval = mb->tag_get_handle( GLOBAL_ID_NAME, 1, MB_TYPE_INTEGER, gid );
  CHECK_ERR(rval);

  Tag conn_ids;
  rval = mb->tag_get_handle( CONN_IDS_NAME, 4, MB_TYPE_INTEGER, conn_ids );
  CHECK_ERR(rval);
  
  std::vector<EntityHandle> quads;
  rval = mb->get_entities_by_type( 0, MBQUAD, quads );
  CHECK_ERR(rval);
  
  std::vector<EntityHandle> conn;
  rval = mb->get_connectivity( &quads[0], quads.size(), conn, true );
  CHECK_ERR(rval);
  
  CHECK_EQUAL( 4*quads.size(), conn.size() );
  std::vector<int> exp_ids(4*quads.size()), act_ids(4*quads.size());
  rval = mb->tag_get_data( conn_ids, &quads[0], quads.size(), &exp_ids[0] );
  CHECK_ERR(rval);
  rval = mb->tag_get_data( gid, &conn[0], conn.size(), &act_ids[0] );
  CHECK_ERR(rval);

  CHECK_EQUAL( exp_ids, act_ids );
}

void check_meshset_common( bool ordered )
{
  ErrorCode rval;
  
  Tag set_ids;
  rval = mb->tag_get_handle( SET_IDS_NAME, ENTS_PER_SET, MB_TYPE_INTEGER, set_ids );
  CHECK_ERR(rval);
  
  Tag gid;
  rval = mb->tag_get_handle( GLOBAL_ID_NAME, 1, MB_TYPE_INTEGER, gid );
  CHECK_ERR(rval);
  
  Range sets;
  rval = mb->get_entities_by_type_and_tag( 0, MBENTITYSET, &set_ids, 0, 1, sets );
  CHECK_ERR(rval);
  CHECK(!sets.empty());
  
  EntityHandle set = 0;
  unsigned flags;
  for (Range::iterator it = sets.begin(); it != sets.end(); ++it) {
    rval = mb->get_meshset_options( *it, flags );
    CHECK_ERR(rval);
    if (( ordered &&  (flags & MESHSET_ORDERED)) ||
        (!ordered && !(flags & MESHSET_ORDERED))) {
      set = *it;
      break;
    }
  }
  CHECK(0 != set);
  
  std::vector<EntityHandle> ents;
  rval = mb->get_entities_by_handle( set, ents );
  CHECK_ERR(rval);
  CHECK_EQUAL( ENTS_PER_SET,(int)ents.size() );
  
  int exp[ENTS_PER_SET], act[ENTS_PER_SET];
  rval = mb->tag_get_data( set_ids, &set, 1, exp );
  CHECK_ERR(rval);
  rval = mb->tag_get_data( gid, &ents[0], ENTS_PER_SET, act );
  CHECK_ERR(rval);
  
  if (!ordered) {
    std::sort( exp, exp+ENTS_PER_SET );
    std::sort( act, act+ENTS_PER_SET );
  }
  
  CHECK_ARRAYS_EQUAL( exp, ENTS_PER_SET, act, ENTS_PER_SET );
  
  if (!(flags & MESHSET_TRACK_OWNER))
    return;
    
  for (int i = 0; i < ENTS_PER_SET; ++i) {
    std::vector<EntityHandle> adj;
    rval = mb->get_adjacencies( &ents[i], 1, 4, false, adj );
    CHECK_ERR(rval);
    CHECK( std::find(adj.begin(), adj.end(), set) != adj.end() );
  }
}

void check_set_meshset()
{
  check_meshset_common(false);
}

void check_list_meshset()
{
  check_meshset_common(true);
}

void check_big_meshset()
{
    // Mesh should have a single set that contains all the vertices.
    // Find it.
  Range sets;
  ErrorCode rval = mb->get_entities_by_type( 0, MBENTITYSET, sets );
  CHECK_ERR(rval);
  
  Range verts;
  rval = mb->get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);  
  
  bool found = false;
  for (Range::iterator it = sets.begin(); it != sets.end(); ++it) {
    Range ents;
    rval = mb->get_entities_by_handle( *it, ents );
    CHECK_ERR(rval);
    if (ents == verts) {
      found = true;
      break;
    }
  }
  CHECK(found);
}


void check_handle_tag()
{
  Range::iterator it;
  ErrorCode rval;

    // first check tag on sets, for which the values should have been
    // updated according to the reordering

  Tag set_handles;
  rval = mb->tag_get_handle( SET_HANDLES_NAME, ENTS_PER_SET, MB_TYPE_HANDLE, set_handles );
  CHECK_ERR(rval);
  
  Range sets;
  rval = mb->get_entities_by_type_and_tag( 0, MBENTITYSET, &set_handles, 0, 1, sets );
  CHECK_ERR(rval);
  CHECK(!sets.empty());

  for (it = sets.begin(); it != sets.end(); ++it) {
    std::vector<EntityHandle> ents;
    rval = mb->get_entities_by_handle( *it, ents );
    CHECK_ERR(rval);
    
    std::vector<EntityHandle> handles(ENTS_PER_SET);
    rval = mb->tag_get_data( set_handles, &*it, 1, &handles[0] );
    CHECK_ERR(rval);
    
    unsigned flags;
    rval = mb->get_meshset_options( *it, flags );
    CHECK_ERR(rval);
    if (!(flags & MESHSET_ORDERED)) 
      std::sort( handles.begin(), handles.end() );
    
    CHECK_EQUAL( ents, handles );
  }
  
    // Now check handle tag on quads.  This tag need to both be re-ordered
    // and have the contained handles updated.
    
  Tag conn_handles;
  rval = mb->tag_get_handle( CONN_NAME, 4, MB_TYPE_HANDLE, conn_handles );
  CHECK_ERR(rval);

  std::vector<EntityHandle> quads;
  rval = mb->get_entities_by_type( 0, MBQUAD, quads );
  CHECK_ERR(rval);
  
  std::vector<EntityHandle> conn;
  rval = mb->get_connectivity( &quads[0], quads.size(), conn, true );
  CHECK_ERR(rval);
  
  std::vector<EntityHandle> tagvals(4*quads.size());
  rval = mb->tag_get_data( conn_handles, &quads[0], quads.size(), &tagvals[0] );
  CHECK_ERR(rval);

  CHECK_EQUAL( conn, tagvals );
}

void check_varlen_tag()
{
  ErrorCode rval;
  
  Tag gid;
  rval = mb->tag_get_handle( GLOBAL_ID_NAME, 1, MB_TYPE_INTEGER, gid );
  CHECK_ERR(rval);
    
  Tag var_data;
  rval = mb->tag_get_handle( VAR_INTS_NAME, 0, MB_TYPE_INTEGER, var_data );
  CHECK_ERR(rval);

  Range verts;
  rval = mb->get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  
  std::vector<int> gids(verts.size());
  rval = mb->tag_get_data( gid, verts, &gids[0] );
  CHECK_ERR(rval);
  
  std::vector<const void*> ptrs(verts.size());
  std::vector<int> sizes(verts.size());
  rval = mb->tag_get_by_ptr( var_data, verts, &ptrs[0], &sizes[0] );
  CHECK_ERR(rval);
  
  for (size_t i = 0; i < gids.size(); ++i) {
    std::vector<int> exp;
    tag_vals_from_gid( gids[i], exp );
    CHECK_ARRAYS_EQUAL( &exp[0], exp.size(), (const int*)ptrs[i], sizes[i] );
  }
}


void check_bit_tag()
{
  ErrorCode rval;
  
  Tag gid;
  rval = mb->tag_get_handle( GLOBAL_ID_NAME, 1, MB_TYPE_INTEGER, gid );
 
  Tag bit_data;
  rval = mb->tag_get_handle( BIT_NAME, BITS_PER_TAG, MB_TYPE_BIT, bit_data );
  CHECK_ERR(rval);

  Range verts;
  rval = mb->get_entities_by_type( 0, MBVERTEX, verts );
  CHECK_ERR(rval);
  
  std::vector<int> gids(verts.size());
  rval = mb->tag_get_data( gid, verts, &gids[0] );
  CHECK_ERR(rval);

  std::vector<unsigned char> exp(gids.size()), act(gids.size());
  for (size_t i = 0; i < exp.size(); ++i)
    exp[i] = bits_from_gid( gids[i] );
  
  rval = mb->tag_get_data( bit_data, verts, &act[0] );
  CHECK_ERR(rval);
  
  CHECK_EQUAL( exp, act );
}
