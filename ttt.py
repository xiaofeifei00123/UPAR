#!/home/fengxiang/anaconda3/envs/wrfout/bin/python
# -*- encoding: utf-8 -*-
'''
Description:

-----------------------------------------
Time             :2021/06/24 18:50:05
Author           :Forxd
Version          :1.0
'''


from typing import Mapping


class A():

    def get_aa(self,):
        # print("yes")
        cc = B()
        cc.get_bb()
        # self.get_bb()
        pass
class B():
    """[summary]
    """
    pass

    def get_bb(self,):
        print("hello")
if __name__ == '__main__':
    # main()
    pass
    aa = A()
    aa.get_aa()
    

