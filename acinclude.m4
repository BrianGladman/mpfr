
AC_DEFUN(AC_MY_LIBS,
[
if ` test "$1" `
then  
  AC_MSG_CHECKING(gmp library)
	if  test -r "$1/lib$2.a"
	then
	  LDADD="$LDADD $1/lib$2.a"
	else
	   AC_MSG_ERROR($2 not found)
	fi
  AC_MSG_RESULT(yes)
else
  AC_CHECK_LIB($2, main, , AC_MSG_ERROR($2 not found))
fi
]
)

AC_DEFUN(AC_MY_HEADERS, 
[
if  test "$1" 
then  
  AC_CHECK_HEADER($1/$2, INCLUDES="$INCLUDES -I$1",AC_MSG_ERROR(echo $2 not found in $1)) 
else
  AC_CHECK_HEADER($2,, 	  AC_MSG_ERROR($2 not found))
fi
])

AC_DEFUN(AC_CHECK_OS, 
[
	AC_MSG_CHECKING(OS type)
	OS_TYPE=`uname -a | awk '{print $ 1}' `
	AC_MSG_RESULT($OS_TYPE)
])

AC_DEFUN(AC_CHECK_MACHTYPE,
[
	AC_MSG_CHECKING(Mach type)
	MACHTYPE=`uname -m`
	AC_MSG_RESULT($MACHTYPE)
])