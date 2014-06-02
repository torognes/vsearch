#include "vsearch.h"


/*
  TODO: Extend it to test query.cc as well, and then do a diff between
  the output of bzquery.cc and query.cc
*/


int main(int argc, char * argv[])
{
  if (argc != 2)
   {
     printf ("syntax: %s [BZ2 FILE]\n", argv[0]);
     return (1);
   }

  char * head;
  char * seq;
  long headlen;
  long seqlen;
  long qno;

  query_bz_open(argv[1]);

  while (query_bz_getnext (&head, &headlen, &seq, &seqlen, &qno))
    printf ("Header: %s   Len: %ld\nSequence: %s   Len: %ld\nId: %ld\n\n", head, headlen, seq, seqlen, qno);

  query_bz_close();

  return (0);
}
