import numpy as np 
import os

def checkMakeDir(path):
    if not os.path.exists(path):
        os.makedirs(path)

class ListDict(object):
    def __init__(self):
        self.item_to_position = {}
        self.items = []

    def __contains__(self, item):
        return item in self.item_to_position

    def __iter__(self):
        return iter(self.items)

    def __len__(self):
        return len(self.items)

    def add(self, item):
        if item in self.item_to_position:
            return
        self.items.append(item)
        self.item_to_position[item] = len(self.items)-1

    def remove(self, item):
        if not item in self.item_to_position:
            return 
        position = self.item_to_position.pop(item)
        last_item = self.items.pop()
        if position != len(self.items):
            self.items[position] = last_item
            self.item_to_position[last_item] = position