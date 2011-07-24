# define m4rie_check(expr)						\
  if (!expr) {								\
    fail_ret += 1;                                                      \
    printf("%s in %s:%d failed\n",__STRING(expr), __FILE__, __LINE__);  \
  } 
