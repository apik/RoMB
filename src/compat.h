#ifndef __COMPAT_H__
#define __COMPAT_H__
#ifndef BOOST_ASSERT_MSG
#define BOOST_ASSERT_MSG(a,b) assert((a)&&(b));
#endif
#endif