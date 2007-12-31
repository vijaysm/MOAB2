#include "MBSimplexTemplateTagAssigner.hpp"

#include "MBMeshRefiner.hpp"

/// Construct a template tag assigner.
MBSimplexTemplateTagAssigner::MBSimplexTemplateTagAssigner( MBMeshRefiner* )
{
}

/// Empty destructor for good form.
MBSimplexTemplateTagAssigner::~MBSimplexTemplateTagAssigner()
{
}

void MBSimplexTemplateTagAssigner::operator()( const void* ta, const void* tb, void* tp )
{
  (void)ta;
  (void)tb;
  (void)tp;
}

void MBSimplexTemplateTagAssigner::operator()( const void* ta, const void* tb, const void* tc, void* tp )
{
  (void)ta;
  (void)tb;
  (void)tc;
  (void)tp;
}
