import physics_objects, contact
from physics_objects import Bounds, PhysicsObject
from pygame import Vector2
import math

# QUAD TREE NODE
# Initialized with an anchor position, size, and capacity.
# TODO: implement maximum depth to prevent infinities
class QuadTreeNode(Bounds):
    def __init__(self, pos:Vector2, size:Vector2, capacity:int=4):
        super().__init__(self, pos, size)
        self.capacity = capacity
        self.divided = False
        self.objects = set()

    # PROPERTIES
    # @property
    # def x(self) -> float:
    #     return self.anchor_pos[0]
    # @property
    # def y(self) -> float:
    #     return self.anchor_pos[1]
    # @property
    # def width(self) -> float:
    #     return self.size[0]
    # @property
    # def height(self) -> float:
    #     return self.size[1]
    
    # METHODS
    # Divides node into 4 child nodes.
    def sub_divide(self) -> bool:
        self.divided = True
        # create 4 child nodes
        size = self.size*0.5
        self.nw = QuadTreeNode(self.pos, size, self.capacity)
        self.ne = QuadTreeNode(self.pos+size.x, size, self.capacity)
        self.sw = QuadTreeNode(self.pos+size.y, size, self.capacity)
        self.se = QuadTreeNode(self.pos+size, size, self.capacity)
        # redistribute current objects into child nodes
        for child in self.get_children():
            for obj in self.objects:
                child.insert(obj)
        self.objects = None

    # Insert an object into this node. Returns 'False' if the object is larger than the bounds of this node. 
    # If this node is full, inserts into child node.
    def insert(self, obj:PhysicsObject) -> bool:
        # check if obj is in bounds of node
        if not self.intersects(obj.bounds): return False
        # if not divided, insert into self
        if not self.divided:
            self.objects.add(obj)
            if len(self.objects) > self.capacity: self.sub_divide()
        # if divided, insert into each each child that intersects
        else:
            for child in [self.ne,self.nw,self.se,self.sw]:
                child.insert(obj)

    # Returns nodes that are within the given boundry range.
    def query(self, bounds:Bounds, found:set = set()) -> set:
        if self.intersects(bounds):
            if not self.divided: 
                found.add(self)
            else:
                for child in self.get_children():
                    child.query(bounds, found)
        return found

    # Returns 'True' if the object is found within this node or its children, else returns 'False'.
    def contains(self, obj) -> bool:
        # check if obj is in self
        if not self.divided:
            return (obj in self.objects)
        # if not, check if obj is in child nodes
        else:
            for child in self.get_children():
                if obj in child: return True
            return False

    # Returns 'True' if this node's boundries intersect the given range.
    def intersects(self, bounds:Bounds):
        # check if boundry intersects range
        return (bounds.pos.x < self.pos.x+self.size.x and bounds.pos.y <self.pos.y+self.size.y 
                and self.pos.x < bounds.pos.x+bounds.size.x and self.pos.y < bounds.pos.y+bounds.size.y)

    # Returns a tuple containing this node's children if this node is a branch. Returns 'False' if this node is not a branch.
    def get_children(self) -> bool|tuple:
        if self.divided:
            return (self.ne,self.nw,self.se,self.sw)
        else: return False


# class QuadTreeBranch(Bounds):
#     def __init__(self, pos, size, ne, nw, se, sw):
#         super().__init__(pos, size)
#         self.children = ()


# class QuadTreeLeaf(Bounds):
#     def __init__(self, pos, size):
#         super().__init__(pos, size)





class QuadTree:
    def __init__(self):
        pass


# SPATIAL HASHING TABLE