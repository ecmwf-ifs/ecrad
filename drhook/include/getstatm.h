
/* statm format: (/proc/self/mem under LINUX)

   in units of pages (4096 bytes)

size       total program size
resident   size of in memory portions
shared     number of the pages that are shared
trs        number of pages that are 'code'
drs        number of pages of data/stack
lrs        number of pages of library
dt         number of dirty pages

*/

struct statm
{
    int size;
    int resident;
    int shared;
    int trs;
    int drs;
    int lrs;
    int dt;
};

extern int getstatm(struct statm *sm);

