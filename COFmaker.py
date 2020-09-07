#!/usr/bin/env python
# coding: utf-8

# In[1]:


# imports
from src.COFmaker import *
from data.fragments import fragments
from data.tests import test_dict


# In[2]:


# simple test case: generate the Porph-Ac2 2D COF
_=COF_4fold(fragments,'porphyrin','Ac2',symmetry='tetragonal_xy',save=True)


# In[3]:


# run phthalocyanine test
test=COF_test()
test.run_test(fragments,test_dict,'phthalocyanine_test_v1')


# In[4]:


# cleanup
# test.cleanup(test_dict,'phthalocyanine_test_v1')


# In[ ]:


# run porphyrin test
test=COF_test()
test.run_test(fragments,test_dict,'porphyrin_test_v1')


# In[ ]:


# cleanup
# test.cleanup(test_dict,'porphyrin_test_v1')

