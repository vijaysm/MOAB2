#include <iMesh.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

struct timeval timeval_diff(struct timeval *lhs,struct timeval *rhs)
{
    const int usec = 1000000;
    struct timeval result;

    if(lhs->tv_usec < rhs->tv_usec)
    {
        int delta = (rhs->tv_usec - lhs->tv_usec) / usec + 1;
        rhs->tv_usec -= delta*usec;
        rhs->tv_sec  += delta;
    }

    if(lhs->tv_usec - rhs->tv_usec > usec)
    {
        int *i = 0;
        *i = 5;
    }

    result.tv_sec  = lhs->tv_sec  - rhs->tv_sec;
    result.tv_usec = lhs->tv_usec - rhs->tv_usec;

    return result;
}

char* format_delta(struct timeval *lhs,struct timeval *rhs)
{
    static char str[512];
    struct timeval delta = timeval_diff(lhs,rhs);
    snprintf(str,sizeof(str),"%u.%06u",(unsigned int)delta.tv_sec,
             (unsigned int)delta.tv_usec);
    return str;
}

int main(int argc,char **argv)
{
    int count;
    struct timeval start,end;
    int err;
    iMesh_Instance mesh=0;
    iBase_EntitySetHandle root;
    char *file;
    int i;

    if(argc != 3)
    {
        fprintf(stderr,"Usage: perf count file\n");
        return 1;
    }

    count = atoi(argv[1]);
    if(count == 0)
        return 1;
    file = argv[2];

    /***** 1 *****/
    gettimeofday(&start,0);
    for(i=0; i<count; i++)
    {
        iMesh_Instance mesh;
        iMesh_newMesh("",&mesh,&err,0);
        iMesh_getRootSet(mesh,&root,&err);
        iMesh_load(mesh,root,file,"",&err,strlen(file),0);
        iMesh_dtor(mesh,&err);
    }
    gettimeofday(&end,0);
    printf("%s\n",format_delta(&end,&start));

    iMesh_newMesh("",&mesh,&err,0);
    iMesh_getRootSet(mesh,&root,&err);
    iMesh_load(mesh,root,file,"",&err,strlen(file),0);

    /***** 2 *****/
    iBase_EntityHandle *entities=0;
    int ent_size=0,ent_alloc;
    iBase_EntityHandle *adj=0;
    int adj_size=0,adj_alloc;
    int *indices=0;
    int ind_size=0,ind_alloc;
    int *offsets=0;
    int off_size=0,off_alloc;

    gettimeofday(&start,0);
    for(i=0; i<count; i++)
    {
        free(entities); entities = 0; ent_size = 0;
        free(adj);      adj = 0;      adj_size = 0;
        free(indices);  indices = 0;  ind_size = 0;
        free(offsets);  offsets = 0;  off_size = 0;
            
        iMesh_getAdjEntIndices(mesh,root,iBase_ALL_TYPES,
                               iMesh_ALL_TOPOLOGIES,iBase_ALL_TYPES,
                               &entities,&ent_size,&ent_alloc,
                               &adj,     &adj_size,&adj_alloc,
                               &indices, &ind_size,&ind_alloc,
                               &offsets, &off_size,&off_alloc,
                               &err);
    }
    gettimeofday(&end,0);
    printf("%s\n",format_delta(&end,&start));

    /***** 3 *****/
    {
        iBase_EntityHandle *adj2=0;
        int adj2_size=0,adj2_alloc;
        int *offsets2=0;
        int off2_size=0,off2_alloc;

        gettimeofday(&start,0);
        for(i=0; i<count; i++)
        {
            adj2 = 0;     adj2_size = 0;
            offsets2 = 0; off2_size = 0;
            iMesh_getEntArrAdj(mesh,entities,ent_size,iBase_ALL_TYPES,
                               &adj2,&adj2_size,&adj2_alloc,
                               &offsets2,&off2_size,&off2_alloc,&err);
            free(adj2);
            free(offsets2);
        }
        gettimeofday(&end,0);
        printf("%s\n",format_delta(&end,&start));
    }

    /***** 4 *****/
    {
        iBase_EntityHandle *adj2=0;
        int adj2_size=0,adj2_alloc;
        int *offsets2=0;
        int off2_size=0,off2_alloc;

        gettimeofday(&start,0);
        for(i=0; i<count; i++)
        {
            adj2 = 0;     adj2_size = 0;
            offsets2 = 0; off2_size = 0;
            iMesh_getEntArr2ndAdj(mesh,entities,ent_size,
                                  iBase_EDGE,iBase_VERTEX,
                                  &adj2,&adj2_size,&adj2_alloc,
                                  &offsets2,&off2_size,&off2_alloc,&err);
            free(adj2);
            free(offsets2);
        }
        gettimeofday(&end,0);
        printf("%s\n",format_delta(&end,&start));
    }

    free(entities);
    free(adj);
    free(indices);
    free(offsets);
    iMesh_dtor(mesh,&err);

    return 0;
}
