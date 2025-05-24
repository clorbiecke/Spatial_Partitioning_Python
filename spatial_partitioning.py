import physics_objects, contact
from pygame import Vector2
import math

# QUAD TREE NODE
# Initialized with an anchor position, size, and capacity.
# 
class QuadTreeNode:
    def __init__(self, anchor_pos:tuple[float,float], size:tuple[float,float], capacity:int=4):
        self.anchor_pos = anchor_pos
        self.size = size
        self.capacity = capacity
        self.divided = False

    # PROPERTIES
    @property
    def x(self) -> float:
        return self.anchor_pos[0]
    @property
    def y(self) -> float:
        return self.anchor_pos[1]
    @property
    def width(self) -> float:
        return self.size[0]
    @property
    def height(self) -> float:
        return self.size[1]
    
    # METHODS
    # Divides node into 4 child nodes.
    def sub_divide(self) -> bool:
        if self.divided: return False
        self.divided = True
        # create 4 child nodes
        self.nw = None
        self.ne = None
        self.sw = None
        self.se = None
        # redistribute current objects into child nodes
    
    # Insert an object into this node. Returns 'False' if the object is larger than the bounds of this node. 
    # If this node is full, inserts into child node.
    def insert(self, obj):
        # check if obj is in bounds of node
        # if not divided, insert into self
        # if divided, insert into each each child that intersects 
        pass

    # Returns nodes that are within the given boundry range.
    def query(self, range:tuple[Vector2,Vector2]):
        # 'range' should be the top-left and bottom-right coordinates (or top-left and size)
        pass

    # Returns 'True' if the object is found within this node or its children, else returns 'False'.
    def contains(self, obj) -> bool:
        # check if obj is in self
        # if not, check if obj is in child nodes
        pass

    # Returns 'True' if this node's boundries intersect the given range.
    def intersects(self, range:tuple[Vector2,Vector2]):
        # check if boundry intersects range
        pass
        
    





class QuadTree:
    def __init__(self):
        pass


# SPATIAL HASHING TABLE