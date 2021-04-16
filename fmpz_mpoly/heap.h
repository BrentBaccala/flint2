
#ifndef FMPZ_MPOLY_HEAP_H
#define FMPZ_MPOLY_HEAP_H

typedef struct fmpz_mpoly_heap_s
{
   fmpz_t exp;
   void * next;
} fmpz_mpoly_heap_s;


/* Heap **********************************************************************/

/* defined in mpoly.h */
/* #define HEAP_LEFT(i) (2*(i)) */
/* #define HEAP_RIGHT(i) (2*(i) + 1) */
/* #define HEAP_PARENT(i) ((i)/2) */

#define FMPZ_POLY_HEAP_ASSIGN(h, c1, c2) \
   do {                                \
      fmpz_set(h.exp, c1);             \
      (h).next = (c2);                 \
   } while (0)

MPOLY_INLINE
void * _fmpz_mpoly_heap_pop(fmpz_mpoly_heap_s * heap, slong * heap_len)
{
   fmpz_t exp;
   slong i, j, s = --(*heap_len);
   mpoly_heap_t * x = (mpoly_heap_t *) heap[1].next;

   i = 1;
   j = 2;

   while (j < s)
   {
      if (!__fmpz_gt(heap[j].exp, heap[j + 1].exp))
         j++;
      heap[i] = heap[j];
      i = j;
      j = HEAP_LEFT(j);
   }

   /* insert last element into heap[i] */
   *exp = *(heap[s].exp);
   j = HEAP_PARENT(i);

   while (i > 1 && __fmpz_gt(exp, heap[j].exp))
   {
      heap[i] = heap[j];
      i = j;
      j = HEAP_PARENT(j);
   }

   heap[i] = heap[s];

   return x; 
}

MPOLY_INLINE
int _fmpz_mpoly_heap_insert(fmpz_mpoly_heap_s * heap, fmpz_t exp, void * x,
       slong * next_loc, slong * heap_len)
{
   slong i = *heap_len, j, n = *heap_len;

   if (i != 1 && fmpz_equal(exp, heap[1].exp))
   {
      ((mpoly_heap_t *) x)->next = (mpoly_heap_t *) heap[1].next;
      heap[1].next = x;

      return 0;
   }

   if (*next_loc < *heap_len)
   {
      if (fmpz_equal(exp, heap[*next_loc].exp))
      {
         ((mpoly_heap_t *) x)->next = (mpoly_heap_t *) heap[*next_loc].next;
         heap[*next_loc].next = x;
         return 0;
      }
   }

   while ((j = HEAP_PARENT(i)) >= 1)
   {
      if (!__fmpz_gt(exp, heap[j].exp))
         break;

      i = j;
   }

   if (j >= 1 && fmpz_equal(exp, heap[j].exp))
   {
      ((mpoly_heap_t *) x)->next = (mpoly_heap_t *) heap[j].next;
      heap[j].next = x;
      *next_loc = j;

      return 0;
   }

   (*heap_len)++;

   while (n > i)
   {
      heap[n] = heap[HEAP_PARENT(n)];
      n = HEAP_PARENT(n);
   }

   FMPZ_POLY_HEAP_ASSIGN(heap[i], exp, x);

   return 1;
}

#endif
