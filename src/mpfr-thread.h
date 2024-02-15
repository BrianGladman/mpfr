/* MPFR internal header related to thread-local variables.

Copyright 2005-2024 Free Software Foundation, Inc.
Contributed by the AriC and Caramba projects, INRIA.

This file is part of the GNU MPFR Library.

The GNU MPFR Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MPFR Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MPFR Library; see the file COPYING.LESSER.  If not, see
https://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#ifndef __MPFR_THREAD_H__
#define __MPFR_THREAD_H__

/* Note: If needed, MPFR_NEED_THREAD_LOCK shall be defined before including
   this header. */

/* Note: Let's define MPFR_THREAD_ATTR even after a #error to make the
   error message more visible (e.g. gcc doesn't immediately stop after
   the #error line and outputs many error messages if MPFR_THREAD_ATTR
   is not defined). But some compilers will just output a message and
   may build MPFR "successfully" (without thread support). */
#ifndef MPFR_THREAD_ATTR
# ifdef MPFR_USE_THREAD_SAFE
#  if defined(_MSC_VER)
#   define MPFR_THREAD_ATTR __declspec( thread )
#  elif defined(MPFR_USE_C11_THREAD_SAFE)
#   define MPFR_THREAD_ATTR _Thread_local
#  else
#   define MPFR_THREAD_ATTR __thread
#  endif
# else
#  define MPFR_THREAD_ATTR
# endif
#endif

/* Note: Several versions of some macros are defined below. In order to
   avoid compilation failure that would appear only with some version,
   for macros that are not related to declarations or initializations,
   such a macro must either be an expression with always the same type
   or always be of the "do ... while(0)" form. */

/**************************************************************************/
/**************************************************************************/
/*                    Start of code for thread locking                    */
/**************************************************************************/

/* If MPFR needs a lock mechanism for thread synchro...
   On 2023-04-12, this is currently used only by the shared cache,
   when it is enabled (see mpfr-impl.h). It is better to allow several
   concurrent readers (as this is likely to occur), but the MPFR code
   does not require that, so that a simple mutex is sufficient. Note
   that much time may be required in a write lock (always exclusive),
   as it involves recomputing a value in a higher precision, but the
   time spent in a read lock will just correspond to the copy of a
   rounded value to the destination. So even if a simple mutex is used,
   the wasted time should not be critical. Moreover, since a mutex is
   simpler to implement, this may also be slightly faster in some cases
   (but no tests have been done yet). */
#ifdef MPFR_NEED_THREAD_LOCK

/* TODO: Prefer the POSIX rwlock method as it allows several readers?
   Alternatively, implement a R/W lock in the C11 mutex version.
   See <https://en.wikipedia.org/wiki/Readers%E2%80%93writer_lock>. */

/**************************************************************************/
/*                 ISO C11 thread-locking version (mutex)                 */
/**************************************************************************/

#if defined (MPFR_HAVE_C11_LOCK)
/* NOTE: This version has not been tested. */

#define MPFR_THREAD_LOCK_METHOD "C11"

#include <threads.h>

#define MPFR_LOCK_DECL(_lock)                           \
  mtx_t _lock;

#define MPFR_LOCK_C(E)                                  \
  do {                                                  \
    if ((E) != thrd_success)                            \
      {                                                 \
        fprintf (stderr, "MPFR lock failure\n");        \
        abort ();                                       \
      }                                                 \
  } while (0)

#define MPFR_LOCK_INIT(_lock)    MPFR_LOCK_C(mtx_init(&(_lock), mtx_plain))
#define MPFR_LOCK_CLEAR(_lock)   do { mtx_destroy(&(_lock)); } while (0)
#define MPFR_LOCK_READ(_lock)    MPFR_LOCK_C(mtx_lock(&(_lock)))
#define MPFR_UNLOCK_READ(_lock)  MPFR_LOCK_C(mtx_unlock(&(_lock)))
#define MPFR_LOCK_WRITE(_lock)   MPFR_LOCK_C(mtx_lock(&(_lock)))
#define MPFR_UNLOCK_WRITE(_lock) MPFR_LOCK_C(mtx_unlock(&(_lock)))

#define MPFR_LOCK_READ2WRITE(_lock) do {} while (0)
#define MPFR_LOCK_WRITE2READ(_lock) do {} while (0)

#define MPFR_ONCE_DECL(_once)                           \
  once_flag _once;

#define MPFR_ONCE_INIT_VALUE ONCE_FLAG_INIT

#define MPFR_ONCE_CALL(_once, _func) do {               \
    call_once(&(_once), (_func));                       \
  } while (0)

#define MPFR_NEED_DEFERRED_INIT 1

/**************************************************************************/
/*             POSIX thread-locking version (pthread rwlock)              */
/**************************************************************************/

/* Note: This uses pthread_rwlock_* functions, which is better than
   pthread_mutex_* as it allows one to have several readers, which
   is likely to occur with the MPFR shared cache. */

#elif defined (HAVE_PTHREAD)

#define MPFR_THREAD_LOCK_METHOD "pthread"

#include <pthread.h>

#define MPFR_LOCK_DECL(_lock)                           \
  pthread_rwlock_t _lock;

#define MPFR_LOCK_C(E)                                  \
  do {                                                  \
    if ((E) != 0)                                       \
      {                                                 \
        fprintf (stderr, "MPFR lock failure\n");        \
        abort ();                                       \
      }                                                 \
  } while (0)

#define MPFR_LOCK_INIT(_lock) MPFR_LOCK_C(pthread_rwlock_init(&(_lock), NULL))

#define MPFR_LOCK_CLEAR(_lock) do {                     \
    pthread_rwlock_destroy(&(_lock));                   \
  } while (0)

#define MPFR_LOCK_READ(_lock)    MPFR_LOCK_C(pthread_rwlock_rdlock(&(_lock)))
#define MPFR_UNLOCK_READ(_lock)  MPFR_LOCK_C(pthread_rwlock_unlock(&(_lock)))
#define MPFR_LOCK_WRITE(_lock)   MPFR_LOCK_C(pthread_rwlock_wrlock(&(_lock)))
#define MPFR_UNLOCK_WRITE(_lock) MPFR_LOCK_C(pthread_rwlock_unlock(&(_lock)))

#define MPFR_LOCK_READ2WRITE(_lock) do {                \
    MPFR_UNLOCK_READ(_lock);                            \
    MPFR_LOCK_WRITE(_lock);                             \
  } while (0)

#define MPFR_LOCK_WRITE2READ(_lock) do {                \
    MPFR_UNLOCK_WRITE(_lock);                           \
    MPFR_LOCK_READ(_lock);                              \
  } while (0)

#define MPFR_ONCE_DECL(_once)                           \
  pthread_once_t _once;

#define MPFR_ONCE_INIT_VALUE PTHREAD_ONCE_INIT

#define MPFR_ONCE_CALL(_once, _func) do {               \
    MPFR_LOCK_C(pthread_once (&(_once), (_func)));      \
  } while (0)

#define MPFR_NEED_DEFERRED_INIT 1

/**************************************************************************/
/*       Thread locking needed, but no available/supported methods        */
/**************************************************************************/

#else

/* TODO: Win32 */
# error "No thread lock / unsupported OS."

#endif

/**************************************************************************/
/*                       Thread locking not needed                        */
/**************************************************************************/

#else  /* MPFR_NEED_THREAD_LOCK */

#define MPFR_LOCK_DECL(_lock)
#define MPFR_LOCK_INIT(_lock)       do {} while (0)
#define MPFR_LOCK_CLEAR(_lock)      do {} while (0)
#define MPFR_LOCK_READ(_lock)       do {} while (0)
#define MPFR_UNLOCK_READ(_lock)     do {} while (0)
#define MPFR_LOCK_WRITE(_lock)      do {} while (0)
#define MPFR_UNLOCK_WRITE(_lock)    do {} while (0)
#define MPFR_LOCK_READ2WRITE(_lock) do {} while (0)
#define MPFR_LOCK_WRITE2READ(_lock) do {} while (0)
#define MPFR_ONCE_INIT_VALUE
#define MPFR_ONCE_DECL(_once)
#define MPFR_ONCE_CALL(_once,_func) do {} while (0)

#endif  /* MPFR_NEED_THREAD_LOCK */

/**************************************************************************/
/*                     End of code for thread locking                     */
/**************************************************************************/
/**************************************************************************/

/**************************************************************************/
/**************************************************************************/
/*                    Start of code for deferred init                     */
/**************************************************************************/

/* If MPFR needs a way to init data before using them
   (MPFR_NEED_DEFERRED_INIT may be set by the above code
   for thread locking)... */
#ifdef MPFR_NEED_DEFERRED_INIT

#ifdef MPFR_HAVE_CONSTRUCTOR_ATTR

/*********************** Use constructor extension ************************/
#define MPFR_DEFERRED_INIT_MASTER_DECL(_id, _init, _clear)      \
  __attribute__ ((constructor)) static void                     \
  mpfr_init_cache_ ## _id (void) {                              \
    _init ;                                                     \
  }                                                             \
  __attribute__ ((destructor)) static void                      \
  mpfr_clean_cache_ ## _id (void) {                             \
    _clear;                                                     \
  }

#define MPFR_DEFERRED_INIT_CALL(_master) do {} while (0)
#define MPFR_DEFERRED_INIT_SLAVE_DECL()
#define MPFR_DEFERRED_INIT_SLAVE_VALUE(_id)

#else /* MPFR_HAVE_CONSTRUCTOR_ATTR */

/**************************** Use once semantic ***************************/
#define MPFR_DEFERRED_INIT_MASTER_DECL(_id, _init, _clear)      \
  static void mpfr_once_ ## _id ## _clear_func (void) {         \
    _clear ;                                                    \
  }                                                             \
  static void mpfr_once_ ## _id ## _init_func (void) {          \
    _init;                                                      \
    atexit(mpfr_once_ ## _id ## _clear_func);                   \
  }

/* MPFR_ONCE_CALL() already uses the "do ... while(0)" form. */
#define MPFR_DEFERRED_INIT_CALL(_master)                        \
  MPFR_ONCE_CALL((_master)->once, (_master)->init_once)

#define MPFR_DEFERRED_INIT_SLAVE_DECL()                         \
  MPFR_ONCE_DECL(once)                                          \
  void (*init_once)(void);

#define MPFR_DEFERRED_INIT_SLAVE_VALUE(_id)                     \
  , MPFR_ONCE_INIT_VALUE, mpfr_once_ ## _id ## _init_func

#endif /* MPFR_HAVE_CONSTRUCTOR_ATTR */

#else /* MPFR_NEED_DEFERRED_INIT */

/* No need */
#define MPFR_DEFERRED_INIT_MASTER_DECL(_id_lock, _init, _clear)
#define MPFR_DEFERRED_INIT_CALL(_master) do {} while (0)
#define MPFR_DEFERRED_INIT_SLAVE_DECL()
#define MPFR_DEFERRED_INIT_SLAVE_VALUE(_id)

#endif /* MPFR_NEED_DEFERRED_INIT */

/**************************************************************************/
/*                     End of code for deferred init                      */
/**************************************************************************/
/**************************************************************************/

#endif
