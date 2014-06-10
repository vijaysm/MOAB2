/** @example VisTags.cpp \n
 * \brief tool for visualizing multi level tags  \n
 * <b>To run</b>: VisTags  <inp_file>  <outfile> -O <read_opts> -t <tags> -l <levels>  -d <dim> \n
 *
 * In this example, it is shown how to create some simple tags for those tags that come from 
 *  climate data, multiple levels.
 *  you can read directly nc data, or *.h5m file that will have the tag with multi levels
 *   output will be a vtk file with dense tags of form tag_name_<level> 
 * the tag name might contain a time index too, like T0 or U0
 * <tag> is a list of tags, separated by commas, no spaces
 * <levels> is a list of levels, separated by commas, no spaces
 *  dimension of entities with the tags will be specified with -d (default 2)
 *
 * an example of use
 *
 * VisTags gcrm_rc.nc  out.vtk -O VARIABLE=u -t u0,u1 -l 0,1,2 -d 2
 * (we knew that it had variable u in the file, that it had 256 levels, that there are 2 time
 *  steps, etc)
 *
 * or
 *  VisTags gcrm_rc.nc  out.vtk  -t u0 -l 0,1,2 -d 2
 *  (it will read all variables, but we need to know that u0 will be created as a tag)
 *
 *  the out.vtk file will contain u0_0, u0_1, as simple dense double tags
 */

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

// Include header for MOAB instance and tag conventions for
#include "moab/Core.hpp" 
#include "MBTagConventions.hpp"
#include "moab/FileOptions.hpp"


int main(int argc, char **argv) {

    // instantiate & load a file 
    moab::Interface *mb = new moab::Core();

    moab::ErrorCode rval;
    if (argc <= 1) 
       return 0;

    int dimension = 2;
    char * file_input = argv[1];
    char * file_output = argv[2];
    char * read_opts = NULL;
    char * tags = NULL; // tags to write, separated by commas; it is the name of the tag
    // in moab, it may have index after reading (T0, T1, etc)
    char * levels = NULL; // levels, separated by commas, no spaces ( like 0,1,19 )
    if (argc>3)
    {
      int index=3;
      while (index<argc)
      {
        if (!strcmp( argv[index], "-O")) // this is for reading options, optional
        {
          read_opts=argv[++index];
        }
        if (!strcmp( argv[index], "-t"))
        {
          tags=argv[++index];
        }
        if (!strcmp( argv[index], "-l"))
        {
          levels=argv[++index];
        }
        if (!strcmp( argv[index], "-d"))
        {
          dimension=atoi(argv[++index]);
        }
        index++;
      }
    }
    std::ostringstream opts;
    opts << ";;TAGS=" << tags << ";LEVELS=" << levels << "\0" ;
    moab::FileOptions fo(opts.str().c_str());

    std::vector<std::string> tagsNames;
    std::vector<int>  levelsArray;
    fo.get_strs_option("TAGS", tagsNames);
    fo.get_ints_option("LEVELS", levelsArray);

    // now create double tags for entities of dimension

    rval = mb->load_file(file_input, 0, read_opts);
    if (rval != moab::MB_SUCCESS) {
      std::cout <<"not loading file\n";
      return 1;
    }
    moab::Range ents;
    rval = mb->get_entities_by_dimension(0, dimension, ents);
    if (rval != moab::MB_SUCCESS) {
      std::cout <<"not getting ents\n";
      return 1;
    }
    for (size_t i=0; i<tagsNames.size(); i++)
    {
      std::string tagName = tagsNames[i];
      moab::Tag tagh;
      rval = mb->tag_get_handle(tagName.c_str(), tagh);

      if (rval != moab::MB_SUCCESS) {
        std::cout <<"not getting tag " << tagName.c_str()<<"\n";
        continue;
      }
      int len=0;
      rval = mb->tag_get_length(tagh, len);
      if (rval != moab::MB_SUCCESS) {
        std::cout <<"not getting tag len" << tagName.c_str()<<"\n";
        continue;
      }
      moab::DataType type;
      rval = mb->tag_get_data_type(tagh, type) ;
      if (rval != moab::MB_SUCCESS) {
        std::cout <<"not getting tag type " << tagName.c_str()<<"\n";
        continue;
      }
      int count;
      void * dataptr;// assume double tags, for simplicity
      rval = mb->tag_iterate( tagh,
          ents.begin(),
          ents.end(),
          count,
          dataptr);
      if (rval != moab::MB_SUCCESS || count != (int)ents.size()) {
        std::cout <<"not getting tag iterate right " << tagName.c_str()<<"\n";
        continue;
      }

      // now create a new tag, with a new name, concatenated, and copy data there , for each level
      for (size_t j=0; j<levelsArray.size();j++)
      {
        int level=levelsArray[j];
        if (level >= len)
        {
          std::cout << "level too big at "<< level << "\n";
          continue;
        }
        std::ostringstream newTagName;
        newTagName << tagName <<"_" << level  ;
        moab::Tag newTagh;
        rval = mb->tag_get_handle(newTagName.str().c_str(), 1, type, newTagh,
            moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
        if (rval != moab::MB_SUCCESS ) {
          std::cout <<"not getting new tag " << newTagName.str() <<"\n";
          continue;
        }
        void * newDataPtr;
        rval = mb->tag_iterate( newTagh,
                            ents.begin(),
                            ents.end(),
                            count,
                            newDataPtr);
        if (rval != moab::MB_SUCCESS  || count !=(int) ents.size()) {
          std::cout <<"not getting new tag iterate" << newTagName.str() <<"\n";
          continue;
        }
        if (type==moab::MB_TYPE_DOUBLE)
        {
          double * ptrD = (double*) newDataPtr;
          double *oldData = (double*)dataptr;
          for (int k=0; k<count; k++, ptrD++)
          {
            *ptrD = oldData[level+count*k];
          }
        }
      }
      mb->tag_delete(tagh);// no need for the tag anymore, write it to the new file
    }

    rval = mb->write_file(file_output);
    if (rval != moab::MB_SUCCESS)
      std::cout <<"can't write file " << file_output << "\n";
    else
      std::cout <<"successfully wrote file " << file_output<< "\n";
    return 0;
} 
