#include "TypeSequenceManager.hpp"
#include "EntitySequence.hpp"
#include "SequenceData.hpp"
#include "MBEntityHandle.h"
#include "TestUtil.hpp"

void test_basic();
void test_lower_bound();
void test_upper_bound();
void test_find();
void test_get_entities();
void test_insert_sequence_merge();
void test_insert_sequence_nomerge();
void test_remove_sequence();
void test_replace_subsequence();
void test_erase();
void test_find_free_handle();
void test_find_free_sequence();
void test_is_free_sequence();
void test_is_free_handle();

/* Construct sequence of vertex handles containing 
   the ID ranges: { [3,7], [100,111], [1001] }
   referencing a SequnceData with the range [1,22000]
 */
void make_basic_sequence( TypeSequenceManager& seq );
void make_basic_sequence( TypeSequenceManager& seq,
                          EntitySequence*& seq3to7,
                          EntitySequence*& seq100to111,
                          EntitySequence*& seq1001 );

/* Compare expected sequence contents to actual contents.
 * Also does some consistency checks.
 */
bool seqman_equal( const MBEntityHandle pair_array[][2],
                     unsigned num_pairs,
                     const TypeSequenceManager& seqman );

/*
 * Insert a sequence into a sequence manager.  Delete passed
 * sequence (and optionally SequenceData) if insertion fails.
 */
MBErrorCode insert_seq( TypeSequenceManager& seqman, 
                        MBEntityHandle start_handle,
                        MBEntityID count,
                        SequenceData* data,
                        bool del_data = false );

int main( )
{
  if (RUN_TEST( test_basic )) {
    printf( "BASIC USE TEST FAILED\nCANNOT TEST FURTHER\n" );
    return 1;
  }

  int error_count = 0;
  error_count += RUN_TEST( test_lower_bound );
  error_count += RUN_TEST( test_upper_bound );
  error_count += RUN_TEST( test_find );
  error_count += RUN_TEST( test_get_entities );
  error_count += RUN_TEST( test_insert_sequence_merge );
  error_count += RUN_TEST( test_insert_sequence_nomerge );
  error_count += RUN_TEST( test_remove_sequence );
  error_count += RUN_TEST( test_replace_subsequence );
  error_count += RUN_TEST( test_erase );
  error_count += RUN_TEST( test_find_free_handle );
  error_count += RUN_TEST( test_find_free_sequence );
  error_count += RUN_TEST( test_is_free_sequence );
  
  if (!error_count) 
    printf( "ALL TESTS PASSED\n");
  else
    printf( "%d TESTS FAILED\n", error_count );
  
  return error_count;
}

class DumSeq : public EntitySequence {
  private:
    int valsPerEnt;
  public: 
    DumSeq( MBEntityHandle start, MBEntityID count, SequenceData* data, int vals_per_ent = 0 )
      : EntitySequence( start, count, data), valsPerEnt(vals_per_ent) {}
    DumSeq( SequenceData* data, int vals_per_ent = 0 )
      : EntitySequence( data->start_handle(), data->size(), data), valsPerEnt(vals_per_ent) {}
    
    DumSeq( DumSeq& split_from, MBEntityHandle here )
      : EntitySequence( split_from, here ), valsPerEnt(split_from.valsPerEnt) {}
    
    virtual ~DumSeq() {}
    
    EntitySequence* split( MBEntityHandle here )
      { return new DumSeq(*this, here); }
    
    SequenceData* create_data_subset( MBEntityHandle a, MBEntityHandle b ) const
      { return data()->subset( a, b, 0 ); }
    
    void get_const_memory_use( unsigned long& a, unsigned long& b) const
      { a = b = 0; }
    unsigned long get_per_entity_memory_use( MBEntityHandle, MBEntityHandle) const
      { return 0; }
    
    int values_per_entity() const
      { return valsPerEnt; }
};
  
void make_basic_sequence( TypeSequenceManager& seqman )
{
  EntitySequence *s, *t, *u;
  make_basic_sequence( seqman, s, t, u );
}


MBErrorCode insert_seq( TypeSequenceManager& seqman, 
                        MBEntityHandle start_handle,
                        MBEntityID count,
                        SequenceData* data,
                        bool del_data )
{
  EntitySequence* seq = new DumSeq( start_handle, count, data );
  MBErrorCode rval = seqman.insert_sequence( seq );
  if (MB_SUCCESS != rval) {
    delete seq;
    if (del_data)
      delete data;
  }
  else {
    CHECK( start_handle >= seq->start_handle() );
    CHECK( start_handle+count-1 <= seq->  end_handle() );
  }
  return rval;
}

void make_basic_sequence( TypeSequenceManager& seqman,
                          EntitySequence*& seq1,
                          EntitySequence*& seq2,
                          EntitySequence*& seq3 )
{
  CHECK( seqman.empty() );
  CHECK_EQUAL( (MBEntityID)0, seqman.get_number_entities() );
  
  SequenceData* data = new SequenceData( 0, 1, 22000 );
  CHECK_ERR( insert_seq( seqman,    3,  5, data ) );
  CHECK_ERR( insert_seq( seqman,  100, 12, data ) );
  CHECK_ERR( insert_seq( seqman, 1001,  1, data ) );

  CHECK( !seqman.empty() );
  CHECK_EQUAL( (MBEntityID)18, seqman.get_number_entities() );
  
  TypeSequenceManager::iterator iter = seqman.begin();
  CHECK( iter != seqman.end() );
  seq1 = *iter;
  CHECK_EQUAL( (MBEntityHandle)3, seq1->start_handle() );
  CHECK_EQUAL( (MBEntityHandle)7, seq1->end_handle() );
  CHECK_EQUAL( data, seq1->data() );
  CHECK_EQUAL( (MBEntityID)5, seq1->size() );
  
  ++iter;
  CHECK( iter != seqman.end() );
  seq2 = *iter;
  CHECK_EQUAL( (MBEntityHandle)100, seq2->start_handle() );
  CHECK_EQUAL( (MBEntityHandle)111, seq2->end_handle() );
  CHECK_EQUAL( data, seq2->data() );
  CHECK_EQUAL( (MBEntityID)12, seq2->size() );
  
  ++iter;
  seq3 = *iter;
  CHECK( iter != seqman.end() );
  CHECK_EQUAL( (MBEntityHandle)1001, seq3->start_handle() );
  CHECK_EQUAL( (MBEntityHandle)1001, seq3->end_handle() );
  CHECK_EQUAL( data, seq3->data() );
  CHECK_EQUAL( (MBEntityID)1, seq3->size() );
  
  ++iter;
  CHECK( iter == seqman.end() );
}

void test_basic( )
{
  TypeSequenceManager seqman;
  make_basic_sequence( seqman );
}

void test_lower_bound()
{
  TypeSequenceManager::const_iterator i;
  TypeSequenceManager seqman;
  EntitySequence *seq1, // 3 to 7
                 *seq2, // 100 to 111
                 *seq3; // 1001
  make_basic_sequence( seqman, seq1, seq2, seq3 );
  
  i = seqman.lower_bound( 2 );
  CHECK_EQUAL( seq1, *i );
  i = seqman.lower_bound( 3 );
  CHECK_EQUAL( seq1, *i );
  i = seqman.lower_bound( 4 );
  CHECK_EQUAL( seq1, *i );
  i = seqman.lower_bound( 7 );
  CHECK_EQUAL( seq1, *i );
  
  i = seqman.lower_bound( 8 );
  CHECK_EQUAL( seq2, *i );
  i = seqman.lower_bound( 99 );
  CHECK_EQUAL( seq2, *i );
  i = seqman.lower_bound( 100 );
  CHECK_EQUAL( seq2, *i );
  i = seqman.lower_bound( 110 );
  CHECK_EQUAL( seq2, *i );
  i = seqman.lower_bound( 111 );
  CHECK_EQUAL( seq2, *i );
  
  i = seqman.lower_bound( 112 );
  CHECK_EQUAL( seq3, *i );
  i = seqman.lower_bound( 1000 );
  CHECK_EQUAL( seq3, *i );
  i = seqman.lower_bound( 1001 );
  CHECK_EQUAL( seq3, *i );
  
  i = seqman.lower_bound( 1002 );
  CHECK( i == seqman.end() );
}


void test_upper_bound()
{
  TypeSequenceManager::const_iterator i;
  TypeSequenceManager seqman;
  EntitySequence *seq1, // 3 to 7
                 *seq2, // 100 to 111
                 *seq3; // 1001
  make_basic_sequence( seqman, seq1, seq2, seq3 );
  
  i = seqman.upper_bound( 2 );
  CHECK_EQUAL( seq1, *i );
  
  i = seqman.upper_bound( 3 );
  CHECK_EQUAL( seq2, *i );
  i = seqman.upper_bound( 4 );
  CHECK_EQUAL( seq2, *i );
  i = seqman.upper_bound( 7 );
  CHECK_EQUAL( seq2, *i );
  i = seqman.upper_bound( 8 );
  CHECK_EQUAL( seq2, *i );
  i = seqman.upper_bound( 99 );
  CHECK_EQUAL( seq2, *i );
  
  i = seqman.upper_bound( 100 );
  CHECK_EQUAL( seq3, *i );
  i = seqman.upper_bound( 110 );
  CHECK_EQUAL( seq3, *i );
  i = seqman.upper_bound( 111 );
  CHECK_EQUAL( seq3, *i );
  i = seqman.upper_bound( 112 );
  CHECK_EQUAL( seq3, *i );
  i = seqman.upper_bound( 1000 );
  CHECK_EQUAL( seq3, *i );

  i = seqman.upper_bound( 1001 );
  CHECK( i == seqman.end() );
}

void test_find()
{
  TypeSequenceManager seqman;
  EntitySequence *seq1, // 3 to 7
                 *seq2, // 100 to 111
                 *seq3, // 1001
                 *seq;
  make_basic_sequence( seqman, seq1, seq2, seq3 );
  
  seq = seqman.find( 2 );
  CHECK_EQUAL( NULL, seq );
  seq = seqman.find( 3 );
  CHECK_EQUAL( seq1, seq );
  seq = seqman.find( 4 );
  CHECK_EQUAL( seq1, seq );
  seq = seqman.find( 7 );
  CHECK_EQUAL( seq1, seq );
  seq = seqman.find( 8 );
  CHECK_EQUAL( NULL, seq );

  seq = seqman.find( 99 );
  CHECK_EQUAL( NULL, seq );
  seq = seqman.find( 100 );
  CHECK_EQUAL( seq2, seq );
  seq = seqman.find( 110 );
  CHECK_EQUAL( seq2, seq );
  seq = seqman.find( 111 );
  CHECK_EQUAL( seq2, seq );
  seq = seqman.find( 112 );
  CHECK_EQUAL( NULL, seq );

  seq = seqman.find( 1000 );
  CHECK_EQUAL( NULL, seq );
  seq = seqman.find( 1001 );
  CHECK_EQUAL( seq3, seq );
  seq = seqman.find( 1002 );
  CHECK_EQUAL( NULL, seq );
}

bool seqman_equal( const MBEntityHandle pair_array[][2],
                     unsigned num_pairs,
                     const TypeSequenceManager& seqman )
{
  unsigned i;
  TypeSequenceManager::const_iterator j = seqman.begin();
  EntitySequence* seq = 0;
  for (i = 0; i < num_pairs; ++i, ++j) {
    if (j == seqman.end())
      break;
    
    if (seq && seq->end_handle() >= (*j)->start_handle()) {
      printf( "Sequence [%lu,%lu] overlaps sequence [%lu,%lu]\n",
       (unsigned long)ID_FROM_HANDLE( seq->start_handle()),
       (unsigned long)ID_FROM_HANDLE( seq->  end_handle()),
       (unsigned long)ID_FROM_HANDLE((*j)->start_handle()),
       (unsigned long)ID_FROM_HANDLE((*j)->  end_handle()));
      return false;
    }
    
    if (seq && seq->data() != (*j)->data() && 
        seq->data()->end_handle() >= (*j)->data()->start_handle()) {
      printf( "SequenceData [%lu,%lu] overlaps SequenceData [%lu,%lu]\n",
       (unsigned long)ID_FROM_HANDLE( seq->data()->start_handle()),
       (unsigned long)ID_FROM_HANDLE( seq->data()->  end_handle()),
       (unsigned long)ID_FROM_HANDLE((*j)->data()->start_handle()),
       (unsigned long)ID_FROM_HANDLE((*j)->data()->  end_handle()));
      return false;
    }
    
    seq = *j;
    if (seq->start_handle() > seq->end_handle()) {
      printf( "Inverted sequence [%lu,%lu]\n",
       (unsigned long)ID_FROM_HANDLE( seq->start_handle()),
       (unsigned long)ID_FROM_HANDLE( seq->  end_handle()));
      return false;
    }
   
    if (pair_array[i][0] != seq->start_handle() ||
        pair_array[i][1] != seq->end_handle())
      break;
      
    if (seq->data()->start_handle() > seq->start_handle() ||
        seq->data()->  end_handle() < seq->  end_handle()) {
      printf( "Sequence [%lu,%lu] has data [%lu,%lu]\n",
       (unsigned long)ID_FROM_HANDLE(seq->start_handle()),
       (unsigned long)ID_FROM_HANDLE(seq->  end_handle()),
       (unsigned long)ID_FROM_HANDLE(seq->data()->start_handle()),
       (unsigned long)ID_FROM_HANDLE(seq->data()->  end_handle()));
      return false;
    }
  }
  
  if (i == num_pairs && j == seqman.end())
    return true;
  
  if (i < num_pairs)
    printf("Sequence Mismatch: Expected: [%lu,%lu], got ",
       (unsigned long)ID_FROM_HANDLE(pair_array[i][0]),
       (unsigned long)ID_FROM_HANDLE(pair_array[i][1]) );
  else
    printf("Sequence Mismatch: Expected END, got " );
  
  if (j == seqman.end())
    printf( "END.\n" );
  else
    printf("[%lu,%lu]\n",
      (unsigned long)ID_FROM_HANDLE((*j)->start_handle()),
      (unsigned long)ID_FROM_HANDLE((*j)->end_handle()) );
  
  return false;
}

void test_get_entities()
{
  TypeSequenceManager seqman;
  make_basic_sequence( seqman );

  CHECK( !seqman.empty() );
  CHECK_EQUAL( (MBEntityID)18, seqman.get_number_entities() );

  MBRange entities;
  seqman.get_entities( entities );
  CHECK_EQUAL( (MBEntityHandle)18, entities.size() );

  MBEntityHandle pairs[][2] = { {3, 7}, {100, 111}, {1001, 1001} };
  CHECK( seqman_equal( pairs, 3, seqman ) );
}


void test_insert_sequence_merge()
{
  TypeSequenceManager seqman;
  make_basic_sequence( seqman );
  SequenceData* data = (*seqman.begin())->data();

    // append a sequence
  CHECK_ERR( insert_seq( seqman, 1003, 1, data ) );
  MBEntityHandle exp1[][2] = { { 3, 7}, {100, 111}, {1001, 1001}, {1003, 1003} };
  CHECK( seqman_equal( exp1, 4, seqman ) );
  
    // prepend a sequence
  CHECK_ERR( insert_seq( seqman, 1, 1, data ) );
  MBEntityHandle exp2[][2] = { {1, 1}, { 3, 7}, {100, 111}, {1001, 1001}, {1003, 1003} };
  CHECK( seqman_equal( exp2, 5, seqman ) );
  
    // insert sequence in middle
  CHECK_ERR( insert_seq( seqman, 150, 11, data ) );
  MBEntityHandle exp3[][2] = { {1, 1}, { 3, 7}, {100, 111}, {150,160}, {1001, 1001}, {1003, 1003} };
  CHECK( seqman_equal( exp3, 6, seqman ) );
  
    // merge sequence with predecessor
  CHECK_ERR( insert_seq( seqman, 8, 13, data ) );
  MBEntityHandle exp4[][2] = { {1, 1}, { 3, 20}, {100, 111}, {150,160}, {1001, 1001}, {1003, 1003} };
  CHECK( seqman_equal( exp4, 6, seqman ) );
  
    // merge sequence with following one
  CHECK_ERR( insert_seq( seqman, 87, 13, data ) );
  MBEntityHandle exp5[][2] = { {1, 1}, { 3, 20}, {87, 111}, {150,160}, {1001, 1001}, {1003, 1003} };
  CHECK( seqman_equal( exp5, 6, seqman ) );
  
    // merge sequence with two adjacent ones
  CHECK_ERR( insert_seq( seqman, 2, 1, data ) );
  MBEntityHandle exp6[][2] = { {1, 20}, {87, 111}, {150,160}, {1001, 1001}, {1003, 1003} };
  CHECK( seqman_equal( exp6, 5, seqman ) );
  
    // try to insert a sequence that overlaps on the end
  CHECK_EQUAL( MB_ALREADY_ALLOCATED, insert_seq( seqman, 900, 1001, data ) );
  
    // try to insert a sequence that overlaps at the start
  CHECK_EQUAL( MB_ALREADY_ALLOCATED, insert_seq( seqman, 111, 140, data ) );
}


void test_insert_sequence_nomerge()
{
  TypeSequenceManager seqman;
  
    // make sure inserting a sequence w/out a SequenceData fails
  CHECK_EQUAL( MB_FAILURE, insert_seq( seqman, 1, 5, NULL ) );
  CHECK( seqman.empty() );
  
    // Now set up a TypeSequenceManager for testing
  
    // Insert an EntitySequence for which the corresponding SequenceData
    // is exactly the same size.
  SequenceData *data1 = new SequenceData( 0, 3, 7 );
  DumSeq* seq = new DumSeq( data1 );
  MBErrorCode rval = seqman.insert_sequence( seq );
  if (MB_SUCCESS != rval) {
    delete seq;
    delete data1;
  }
  CHECK_ERR( rval );
  
    // Insert an EntitySequence with additional room on both ends of
    // the SequenceData
  SequenceData *data2 = new SequenceData( 0, 100, 999 );
  seq = new DumSeq( 200, 100, data2 );
  rval = seqman.insert_sequence( seq );
  if (MB_SUCCESS != rval) {
    delete seq;
    delete data2;
  }
  CHECK_ERR( rval );
  
    // Insert another EntitySequence sharing the previous SequenceData
  seq = new DumSeq( 400, 100, data2 );
  rval = seqman.insert_sequence( seq );
  if (MB_SUCCESS != rval)
    delete seq;
  CHECK_ERR( rval );
  
    
      // Setup complete, begin tests
  
    // Test inserting sequence that appends an existing sequence
    // but overlaps underling SequenceData boundary
  SequenceData* data = new SequenceData( 0, 999, 1000 );
  seq = new DumSeq( data );
  rval = seqman.insert_sequence( seq );
  CHECK_EQUAL( MB_ALREADY_ALLOCATED, rval );
  delete seq;
  delete data;
  
    // Test inserting sequence that prepends an existing sequence
    // but overlaps underling SequenceData boundary
  data = new SequenceData( 0, 50, 199 );
  seq = new DumSeq( data );
  rval = seqman.insert_sequence( seq );
  CHECK_EQUAL( MB_ALREADY_ALLOCATED, rval );
  delete seq;
  delete data;
  
    // Test fits within existing, but has different SequenceData
  data = new SequenceData( 0, 500, 599 );
  seq = new DumSeq( data );
  rval = seqman.insert_sequence( seq );
  CHECK_EQUAL( MB_ALREADY_ALLOCATED, rval );
  delete seq;
  delete data;
  
    // Make sure we're starting out with what we exoect
  MBEntityHandle exp1[][2] = { { 3, 7}, {200, 299}, {400, 499} };
  CHECK( seqman_equal( exp1, 3, seqman ) );
  
    // Test fits within existing, and has same data
  CHECK_ERR( insert_seq( seqman, 600, 100, data2 ) );
  MBEntityHandle exp2[][2] = { { 3, 7}, {200, 299}, {400, 499}, {600, 699} };
  CHECK( seqman_equal( exp2, 4, seqman ) );
  
    // Test is entirely outside existing data
  CHECK_ERR( insert_seq( seqman, 2000, 2, new SequenceData( 0, 2000, 2001 ), true ) );
  MBEntityHandle exp3[][2] = { { 3, 7}, {200, 299}, {400, 499}, {600, 699}, {2000, 2001} };
  CHECK( seqman_equal( exp3, 5, seqman ) );
  
    // Test abutts end of exising data
  CHECK_ERR( insert_seq( seqman, 1000, 6, new SequenceData( 0, 1000, 1005 ), true ) );
  MBEntityHandle exp4[][2] = { { 3, 7}, {200, 299}, {400, 499}, {600, 699}, {1000,1005}, {2000, 2001} };
  CHECK( seqman_equal( exp4, 6, seqman ) );
 
    // Test abutts beginning of exising data
  CHECK_ERR( insert_seq( seqman, 50, 50, new SequenceData( 0, 50, 99 ), true ) );
  MBEntityHandle exp5[][2] = { { 3, 7}, {50, 99}, {200, 299}, {400, 499}, {600, 699}, {1000,1005}, {2000, 2001} };
  CHECK( seqman_equal( exp5, 7, seqman ) );
}

void test_remove_sequence()
{
  TypeSequenceManager seqman;
  EntitySequence *seq1, // 3 to 7
                 *seq2, // 100 to 111
                 *seq3; // 1001
  make_basic_sequence( seqman, seq1, seq2, seq3 );

    // test removing something that hasn't been inserted
  bool last;
  DumSeq junk( 3, 5, NULL );
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, seqman.remove_sequence( &junk, last ) );
  MBEntityHandle exp1[][2] = { { 3, 7}, {100, 111}, {1001, 1001} };
  CHECK( seqman_equal( exp1, 3, seqman ) );
  
    // remove the middle sequence
  CHECK_ERR( seqman.remove_sequence( seq2, last ) );
  CHECK( !last );
  delete seq2;
  
    // remove the first sequence
  CHECK_ERR( seqman.remove_sequence( seq1, last ) );
  CHECK( !last );
  delete seq1;
  
    // remove the last sequence
  CHECK_ERR( seqman.remove_sequence( seq3, last ) );
  CHECK( last );
  SequenceData* data = seq3->data();
  delete seq3;
  delete data;
}

void test_replace_subsequence()
{
  MBErrorCode rval;
  TypeSequenceManager seqman;

    // create an intial set
  SequenceData* data1 = new SequenceData( 0, 51, 950 );
  CHECK_ERR( insert_seq( seqman, 101, 100, data1, true ) );
  CHECK_ERR( insert_seq( seqman, 301, 300, data1 ) );
  CHECK_ERR( insert_seq( seqman, 701, 100, data1 ) );
  
    // try a sequence that is outside all existing data
  SequenceData* data = new SequenceData( 0, 10, 20 );
  DumSeq* seq = new DumSeq( data );
  rval = seqman.replace_subsequence( seq, 0 );
  CHECK_EQUAL( MB_FAILURE, rval );
  delete seq;
  delete data;
  
    // try a sequence that overlaps the start of the data
  data = new SequenceData( 0, 40, 60 );
  seq = new DumSeq( data );
  rval = seqman.replace_subsequence( seq, 0 );
  CHECK_EQUAL( MB_FAILURE, rval );
  delete seq;
  delete data;
  
    // try a sequence that is within the data but not within any sequence
  data = new SequenceData( 0, 60, 70 );
  seq = new DumSeq( data );
  rval = seqman.replace_subsequence( seq, 0 );
  CHECK_EQUAL( MB_FAILURE, rval );
  delete seq;
  delete data;
  
    // try a sequence that overlaps an existing sequence
  data = new SequenceData( 0, 60, 101 );
  seq = new DumSeq( data );
  rval = seqman.replace_subsequence( seq, 0 );
  CHECK_EQUAL( MB_FAILURE, rval );
  delete seq;
  delete data;
  
    // try a sequence that should work, but with a SequenceData that 
    // overlaps an existing sequence
  data = new SequenceData( 0, 150, 200 );
  seq = new DumSeq( 190, 200, data );
  rval = seqman.replace_subsequence( seq, 0 );
  CHECK_EQUAL( MB_FAILURE, rval );
  delete seq;
  delete data;
  
    // check that we're starting with what we expect
  MBEntityHandle exp1[][2] = { {101, 200}, {301, 600}, {701, 800} };
  CHECK( seqman_equal( exp1, 3, seqman ) );
  
    // split at start of sequence
  data = new SequenceData( 0, 101, 105 );
  seq = new DumSeq( data );
  rval = seqman.replace_subsequence( seq, 0 );
  if (MB_SUCCESS != rval) {
    delete seq;
    delete data;
  }
  CHECK_ERR( rval );
  MBEntityHandle exp2[][2] = { {101, 105}, {106, 200}, {301, 600}, {701, 800} };
  CHECK( seqman_equal( exp2, 4, seqman ) );
  
    // split at end of sequence
  data = new SequenceData( 0, 750, 800 );
  seq = new DumSeq( data );
  rval = seqman.replace_subsequence( seq, 0 );
  if (MB_SUCCESS != rval) {
    delete seq;
    delete data;
  }
  CHECK_ERR( rval );
  MBEntityHandle exp3[][2] = { {101, 105}, {106, 200}, {301, 600}, {701, 749}, {750,800} };
  CHECK( seqman_equal( exp3, 5, seqman ) );
  
    // split at middle of sequence
  data = new SequenceData( 0, 400, 499 );
  seq = new DumSeq( data );
  rval = seqman.replace_subsequence( seq, 0 );
  if (MB_SUCCESS != rval) {
    delete seq;
    delete data;
  }
  CHECK_ERR( rval );
  MBEntityHandle exp4[][2] = { {101, 105}, {106, 200}, {301, 399}, {400,499}, {500,600}, {701, 749}, {750,800} };
  CHECK( seqman_equal( exp4, 7, seqman ) );
}

void test_erase()
{
  TypeSequenceManager seqman;
  make_basic_sequence( seqman );

    // verify initial state
  MBEntityHandle exp1[][2] = { {3, 7}, {100, 111}, {1001, 1001} };
  CHECK( seqman_equal( exp1, 3, seqman ) );

    // try erasing invalid handles at start of existing sequence
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, seqman.erase( 1000, 1001 ) );
    // try erasing invalid entities at end of existing sequence
  CHECK_EQUAL( MB_ENTITY_NOT_FOUND, seqman.erase( 3, 8 ) );
    // verify initial state
  CHECK( seqman_equal( exp1, 3, seqman ) );
  
    // erase from front of sequence
  CHECK_ERR( seqman.erase( 3, 6 ) );
  MBEntityHandle exp2[][2] = { {7, 7}, {100, 111}, {1001, 1001} };
  CHECK( seqman_equal( exp2, 3, seqman ) );
  
    // erase from end of sequence
  CHECK_ERR( seqman.erase( 110, 111 ) );
  MBEntityHandle exp3[][2] = { {7, 7}, {100, 109}, {1001, 1001} };
  CHECK( seqman_equal( exp3, 3, seqman ) );
  
    // erase from middle of sequence
  CHECK_ERR( seqman.erase( 105, 107 ) );
  MBEntityHandle exp4[][2] = { {7, 7}, {100, 104}, {108,109}, {1001, 1001} };
  CHECK( seqman_equal( exp4, 4, seqman ) );
  
    // erase sequence
  CHECK_ERR( seqman.erase( 7, 7 ) );
  MBEntityHandle exp5[][2] = { {100, 104}, {108,109}, {1001, 1001} };
  CHECK( seqman_equal( exp5, 3, seqman ) );
  
    // erase sequence
  CHECK_ERR( seqman.erase( 108, 109 ) );
  MBEntityHandle exp6[][2] = { {100, 104}, {1001, 1001} };
  CHECK( seqman_equal( exp6, 2, seqman ) );
  
    // erase sequence
  CHECK_ERR( seqman.erase( 100, 104 ) );
  MBEntityHandle exp7[][2] = { {1001, 1001} };
  CHECK( seqman_equal( exp7, 1, seqman ) );
  
    // erase sequence
  CHECK_ERR( seqman.erase( 1001, 1001 ) );
  CHECK( seqman.empty() );
}


void test_find_free_handle()
{
  bool append;
  TypeSequenceManager::iterator seq;
  TypeSequenceManager seqman;
  make_basic_sequence( seqman ); // { [3,7], [100,111], [1001] }
  
  seq = seqman.find_free_handle( 0, MB_END_ID, append );
  CHECK( seq != seqman.end() );
    // expect the first available handle (2).
  CHECK_EQUAL( (MBEntityHandle)3, (*seq)->start_handle() );
  CHECK_EQUAL( (MBEntityHandle)7, (*seq)->  end_handle() );
  CHECK( !append );

    // Expect end() if no adjacent sequence
  seq = seqman.find_free_handle( 9, 98, append );
  CHECK( seq == seqman.end() );
  
    // Try a limited handle range
  seq = seqman.find_free_handle( 8, 99, append );
  CHECK( seq != seqman.end() );
    // expect the first available handle (8).
  CHECK_EQUAL( (MBEntityHandle)3, (*seq)->start_handle() );
  CHECK_EQUAL( (MBEntityHandle)7, (*seq)->  end_handle() );
  CHECK( append );
  
    // Try an unambigious case (above tests have multiple
    // possible answers, were we assume the first available
    // handle).
  seq = seqman.find_free_handle( 8, 98, append );
  CHECK( seq != seqman.end() );
  CHECK_EQUAL( (MBEntityHandle)3, (*seq)->start_handle() );
  CHECK_EQUAL( (MBEntityHandle)7, (*seq)->  end_handle() );
  CHECK( append );
  
    // Try an unambigious case (above tests have multiple
    // possible answers, were we assume the first available
    // handle).
  seq = seqman.find_free_handle( 9, 99, append );
  CHECK( seq != seqman.end() );
  CHECK_EQUAL( (MBEntityHandle)100, (*seq)->start_handle() );
  CHECK_EQUAL( (MBEntityHandle)111, (*seq)->  end_handle() );
  CHECK( !append );
  
    // Try a case where the expected result handle
    // is in the middle of the input range.
  seq = seqman.find_free_handle( 900, 1100, append );
  CHECK( seq != seqman.end() );
  CHECK_EQUAL( (MBEntityHandle)1001, (*seq)->start_handle() );
  CHECK_EQUAL( (MBEntityHandle)1001, (*seq)->  end_handle() );
    // Expect first available handle
  CHECK( !append );
}


void test_find_free_sequence()
{
  MBEntityHandle start;
  SequenceData* data = 0;
  MBEntityID data_size = 0;
  TypeSequenceManager seqman;
  make_basic_sequence( seqman ); // { [3,7], [100,111], [1001] }
  SequenceData* expdata = (*seqman.begin())->data();
  
  start = seqman.find_free_sequence( 2, 1, 3, data, data_size );
  CHECK_EQUAL( expdata, data );
  CHECK_EQUAL( (MBEntityHandle)1, start );
  
  start = seqman.find_free_sequence( 3, 1, 7, data, data_size );
  CHECK_EQUAL( (MBEntityHandle)0, start );
  CHECK_EQUAL( NULL, data );
  
  start = seqman.find_free_sequence( 30, 1, 120, data, data_size );
  CHECK_EQUAL( expdata, data );
  CHECK( start == 8 || start == 70 );
  
  start = seqman.find_free_sequence( 10, 92, 999, data, data_size );
  CHECK_EQUAL( expdata, data );
  CHECK_EQUAL( (MBEntityHandle)112, start );
  
  start = seqman.find_free_sequence( 100, 1, 600, data, data_size );
  CHECK_EQUAL( expdata, data );
  CHECK_EQUAL( (MBEntityHandle)112, start );
  
  start = seqman.find_free_sequence( 1000, 1, MB_END_ID, data, data_size );
  CHECK_EQUAL( expdata, data );
  CHECK_EQUAL( (MBEntityHandle)1002, start );
  
  start = seqman.find_free_sequence( 980, 1, 1800, data, data_size );
  CHECK_EQUAL( (MBEntityHandle)0, start );
  CHECK_EQUAL( NULL, data );
}


void test_is_free_sequence()
{
  SequenceData* data = 0;
  TypeSequenceManager seqman;
  make_basic_sequence( seqman ); // { [3,7], [100,111], [1001] }
  SequenceData* expdata = (*seqman.begin())->data();
  
  CHECK( !seqman.is_free_sequence( 1, 3, data ) );
  CHECK(  seqman.is_free_sequence( 1, 2, data ) );
  CHECK_EQUAL( expdata, data );
  CHECK( !seqman.is_free_sequence( 7, 93, data ) );
  CHECK( !seqman.is_free_sequence( 8, 93, data ) );
  CHECK(  seqman.is_free_sequence( 8, 92, data ) );
  CHECK_EQUAL( expdata, data );
  CHECK( !seqman.is_free_sequence( 111, 890, data ) );
  CHECK( !seqman.is_free_sequence( 112, 890, data ) );
  CHECK(  seqman.is_free_sequence( 112, 879, data ) );
  CHECK_EQUAL( expdata, data );
  CHECK(  seqman.is_free_sequence( 1002, 1, data ) );
  CHECK_EQUAL( expdata, data );
  CHECK(  seqman.is_free_sequence( 2000, 20000, data ) );
  CHECK_EQUAL( expdata, data );
}    


void test_is_free_handle()
{
  // Construct a TypeSequenceManager with the following data:
  // EntitySequence: |[1,500]|   |[601,1000]|       |[2500,2599]|    |[2800,2999]|
  // SequenceData:   |       [1,1000]       |    |         [2001,3000]            |
  TypeSequenceManager seqman;
  SequenceData* data1 = new SequenceData(0,   1,1000);
  SequenceData* data2 = new SequenceData(0,2001,3000);
  CHECK( insert_seq( seqman,    1, 500, data1, true  ) );
  CHECK( insert_seq( seqman,  601, 400, data1, false ) );
  CHECK( insert_seq( seqman, 2500, 100, data2, true  ) );
  CHECK( insert_seq( seqman, 2800, 200, data2, false ) );
  
  // Begin tests
  TypeSequenceManager::iterator seq;
  SequenceData* data;
  MBEntityHandle first, last;
 
  // Test handle in use

  CHECK_EQUAL( MB_ALREADY_ALLOCATED, seqman.is_free_handle(   1, seq, data, first, last ) );
  CHECK_EQUAL( MB_ALREADY_ALLOCATED, seqman.is_free_handle( 300, seq, data, first, last ) );
  CHECK_EQUAL( MB_ALREADY_ALLOCATED, seqman.is_free_handle( 500, seq, data, first, last ) );
  CHECK_EQUAL( MB_ALREADY_ALLOCATED, seqman.is_free_handle( 601, seq, data, first, last ) );
  CHECK_EQUAL( MB_ALREADY_ALLOCATED, seqman.is_free_handle(2500, seq, data, first, last ) );
  CHECK_EQUAL( MB_ALREADY_ALLOCATED, seqman.is_free_handle(2599, seq, data, first, last ) );
  CHECK_EQUAL( MB_ALREADY_ALLOCATED, seqman.is_free_handle(2800, seq, data, first, last ) );
  CHECK_EQUAL( MB_ALREADY_ALLOCATED, seqman.is_free_handle(2999, seq, data, first, last ) );
 
 
  // Test prepend to sequence
  
  seq = seqman.end(); data = 0; first = last = 0;
  CHECK( seqman.is_free_handle( 600, seq, data, first, last ) );
  CHECK( seq != seqman.end() );
  CHECK_EQUAL( (MBEntityHandle)601, (*seq)->start_handle() );
  CHECK_EQUAL( data1, data );
  CHECK_EQUAL( (MBEntityHandle)600, first );
  CHECK_EQUAL( (MBEntityHandle)600, last );
  
  seq = seqman.end(); data = 0; first = last = 0;
  CHECK( seqman.is_free_handle( 2499, seq, data, first, last ) );
  CHECK( seq != seqman.end() );
  CHECK_EQUAL( (MBEntityHandle)2500, (*seq)->start_handle() );
  CHECK_EQUAL( data2, data );
  CHECK_EQUAL( (MBEntityHandle)2499, first );
  CHECK_EQUAL( (MBEntityHandle)2499, last );

  seq = seqman.end(); data = 0; first = last = 0;
  CHECK( seqman.is_free_handle( 2799, seq, data, first, last ) );
  CHECK( seq != seqman.end() );
  CHECK_EQUAL( (MBEntityHandle)2800, (*seq)->start_handle() );
  CHECK_EQUAL( data2, data );
  CHECK_EQUAL( (MBEntityHandle)2799, first );
  CHECK_EQUAL( (MBEntityHandle)2799, last );
 
  // Test append to sequence
  
  seq = seqman.end(); data = 0; first = last = 0;
  CHECK( seqman.is_free_handle( 501, seq, data, first, last ) );
  CHECK( seq != seqman.end() );
  CHECK_EQUAL( (MBEntityHandle)1, (*seq)->start_handle() );
  CHECK_EQUAL( data1, data );
  CHECK_EQUAL( (MBEntityHandle)501, first );
  CHECK_EQUAL( (MBEntityHandle)501, last );
  
  seq = seqman.end(); data = 0; first = last = 0;
  CHECK( seqman.is_free_handle( 2600, seq, data, first, last ) );
  CHECK( seq != seqman.end() );
  CHECK_EQUAL( (MBEntityHandle)2500, (*seq)->start_handle() );
  CHECK_EQUAL( data2, data );
  CHECK_EQUAL( (MBEntityHandle)2600, first );
  CHECK_EQUAL( (MBEntityHandle)2600, last );

  seq = seqman.end(); data = 0; first = last = 0;
  CHECK( seqman.is_free_handle( 3000, seq, data, first, last ) );
  CHECK( seq != seqman.end() );
  CHECK_EQUAL( (MBEntityHandle)2800, (*seq)->start_handle() );
  CHECK_EQUAL( data2, data );
  CHECK_EQUAL( (MBEntityHandle)3000, first );
  CHECK_EQUAL( (MBEntityHandle)3000, last );
  
    // Test new sequence in existing SequenceData

  seq = seqman.end(); data = 0; first = last = 0;
  CHECK( seqman.is_free_handle( 502, seq, data, first, last ) );
  CHECK( seqman.end() == seq );
  CHECK_EQUAL( data1, data );
  CHECK_EQUAL( (MBEntityHandle)501, first );
  CHECK_EQUAL( (MBEntityHandle)600, last );

  seq = seqman.end(); data = 0; first = last = 0;
  CHECK( seqman.is_free_handle( 599, seq, data, first, last ) );
  CHECK( seqman.end() == seq );
  CHECK_EQUAL( data1, data );
  CHECK_EQUAL( (MBEntityHandle)501, first );
  CHECK_EQUAL( (MBEntityHandle)600, last );

  seq = seqman.end(); data = 0; first = last = 0;
  CHECK( seqman.is_free_handle( 2001, seq, data, first, last ) );
  CHECK( seqman.end() == seq );
  CHECK_EQUAL( data2, data );
  CHECK_EQUAL( (MBEntityHandle)2001, first );
  CHECK_EQUAL( (MBEntityHandle)2499, last );

  seq = seqman.end(); data = 0; first = last = 0;
  CHECK( seqman.is_free_handle( 2498, seq, data, first, last ) );
  CHECK( seqman.end() == seq );
  CHECK_EQUAL( data2, data );
  CHECK_EQUAL( (MBEntityHandle)2001, first );
  CHECK_EQUAL( (MBEntityHandle)2499, last );

    // Test new SequenceData 

  seq = seqman.end(); data = 0; first = last = 0;
  CHECK( seqman.is_free_handle( 1001, seq, data, first, last ) );
  CHECK( seqman.end() == seq );
  CHECK_EQUAL( (SequenceData*)0, data );
  CHECK_EQUAL( (MBEntityHandle)1001, first );
  CHECK_EQUAL( (MBEntityHandle)2000, last );

  seq = seqman.end(); data = 0; first = last = 0;
  CHECK( seqman.is_free_handle( 1500, seq, data, first, last ) );
  CHECK( seqman.end() == seq );
  CHECK_EQUAL( (SequenceData*)0, data );
  CHECK_EQUAL( (MBEntityHandle)1001, first );
  CHECK_EQUAL( (MBEntityHandle)2000, last );

  seq = seqman.end(); data = 0; first = last = 0;
  CHECK( seqman.is_free_handle( 2000, seq, data, first, last ) );
  CHECK( seqman.end() == seq );
  CHECK_EQUAL( (SequenceData*)0, data );
  CHECK_EQUAL( (MBEntityHandle)1001, first );
  CHECK_EQUAL( (MBEntityHandle)2000, last );

  seq = seqman.end(); data = 0; first = last = 0;
  CHECK( seqman.is_free_handle( 3001, seq, data, first, last ) );
  CHECK( seqman.end() == seq );
  CHECK_EQUAL( (SequenceData*)0, data );
  CHECK_EQUAL( (MBEntityHandle)3001, first );
  CHECK_EQUAL( (MBEntityHandle)MB_END_ID, last );

  seq = seqman.end(); data = 0; first = last = 0;
  CHECK( seqman.is_free_handle( 10000, seq, data, first, last ) );
  CHECK( seqman.end() == seq );
  CHECK_EQUAL( (SequenceData*)0, data );
  CHECK_EQUAL( (MBEntityHandle)3001, first );
  CHECK_EQUAL( (MBEntityHandle)MB_END_ID, last );
}
