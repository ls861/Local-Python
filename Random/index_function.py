#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 12:59:57 2019

@author: lester
"""

def index(value, list_of_items):
    ''' function that takes a value, finds the closet match in a given list of
    items and then returns the index of where this value is in the list '''
    
    diff = 9999
    index = -1
    
    for i in range(len(list_of_items)):
        if abs(list_of_items[i] - value) < diff:
            index = i
            diff = abs(list_of_items[i] - value)
            
    return index
    
    
    
    
    
    
    
    
value = 5
list_of_items = [8, 40, 11, 0, 3, 6, 3, 8, 9, 2, 3]

print(index(value, list_of_items))


