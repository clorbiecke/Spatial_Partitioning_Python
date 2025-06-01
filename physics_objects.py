import pygame
from pygame import Color
import pygame.gfxdraw
from pygame import Surface

from pygame.math import Vector2, Vector3
# from datetime import datetime
from collections import deque
import colors
import itertools
import random
import math
import copy

class Bounds:
    def __init__(self, pos:Vector2, size:Vector2):
        self.pos = pos
        self.size = size

    @property
    def pos(self) -> Vector2:
        return self._pos
    @pos.setter
    def pos(self, pos:Vector2):
        self._pos = Vector2(pos)

    @property
    def size(self) -> Vector2:
        return self._size
    @size.setter
    def size(self, size:Vector2):
        self._size = Vector2(size)

    @property
    def width(self) -> float:
        return self.size.x
    @property
    def height(self) -> float:
        return self.size.y

class PhysicsObject:
    phys_objs = []
    DRAG_C = 0
    def __init__(self, name="", mass=1, pos=(0,0), vel=(0,0), momi=math.inf, angle=0, avel=0, static:bool=False, generate_contacts=True, on_update:lambda obj:() = lambda obj:()):
        self.name:str = name
        self.static:bool = static
        self.mass:float = mass
        self.pos:Vector2 = pos
        self.vel:Vector2 = Vector2(vel)
        self.force:Vector2 = Vector2(0,0)
        self.momi:float = momi
        self.angle:float = angle # in degrees
        self.avel:float = avel # deg/sec
        self.torque:float = 0 # rotational force
        self.set_bounds()
    
        self.generate_contacts:bool = generate_contacts

        self.phys_objs.append(self)
        self.on_update = on_update
        self.on_collision_enter = []
        self.colliding_objects = []
    
    def __del__(self):
        # print(f"{self.name} is being deleted at {datetime.now()}...")
        pass

    #-- METHODS
    # Sets the bounds of this object. Should be overridden in subclasses.
    def set_bounds(self):
        self.bounds:Bounds = Bounds(self.pos.x, self.pos.y, 0,0)
    # Updates the bounds of this object. Should be overridden in subclasses.
    def update_bounds(self):
        self.bounds.pos = self.pos
    # Updates the shape of this object, then calls update_bounds(). Called at the end of update() before on_update(). Should be overridden in subclasses.
    def update_shape(self, dt=0):
        self.update_bounds()

    def clear_force(self):
        self.force = Vector2(0,0)
        self.torque = 0

    def add_force(self, force):
        self.force += Vector2(force)

    def set(self, pos=None, angle=None):
        if pos is not None:
            self.pos = Vector2(pos)
        if angle is not None:
            self.angle = angle
    
    def add_torque(self, torque):
        self.torque += torque
    
    # Adds impulse, updating velocity. (impulse = force * delta_time)
    def add_impulse(self, impulse:Vector2, point:Vector2=None):
        # Delta v = J / m
        if point is not None:
            S:Vector2 = point - self.pos
            # delta avel = (S x J) / I
            #NOTE: is this right? used to be avel = _ rather than avel += _
            d_avel = math.degrees(S.cross(impulse) / self.momi)
            self.avel += d_avel
        self.vel += impulse / self.mass

    # Change position by adding delta.
    def delta_pos(self, delta:Vector2):
        # change position by delta
        self.pos += delta

    # call to update the moment of inertia if mass or center of mass changes
    def update_momi(self):
        pass

    def update(self, dt=0):
        # only update if not static
        if not self.static:
            # update velocity using the current force
            self.vel += (self.force/self.mass) * dt
            # update position using the newly updated velocity
            self.pos += (self.vel * dt)
            # round vel to 0 if its low enough
            if (self.vel.magnitude() < 0.001): self.vel = Vector2(0)
            # update angular velocity
            self.avel += (self.torque/self.momi) * dt
            # update angle
            self.angle += (self.avel * dt)
            # update shape
            self.update_shape(dt)

        
        for fun in self.on_update:
            fun(dt)

    def collision_enter(self, other, **kwargs):
        kwargs = kwargs
        for fun in self.on_collision_enter:
            fun(other)


    # PROPERTIES
    @property
    def name(self):
        return self._name() if callable(self._name) else self._name
    @name.setter
    def name(self, name:str):
        if name == "":
            name = self.__class__.__name__
        self._name = name
    
    @property
    def mass(self):
        if self.static: 
            return math.inf
        else:
            return float(self._mass()) if callable(self._mass) else float(self._mass)
    @mass.setter
    def mass(self, mass:float):
        self._mass = mass
        self.update_momi()

    @property
    def area(self):
        return 0

    @property
    def pos(self) -> Vector2:
        return self._pos() if callable(self._pos) else self._pos
    @pos.setter
    def pos(self, pos:Vector2):
        self._pos = Vector2(pos)

    @property
    def on_update(self) -> list[lambda:()]:
        return self._on_update
    @on_update.setter
    def on_update(self, funs:list[lambda:()]):
        if isinstance(funs, list):
            self._on_update = funs
        else: self._on_update = [funs]

    # TODO: add object to list? make sure this is actually working (add call to this in contact.py)
    @property
    def on_collision_enter(self) -> list[lambda:()]:
        return self._on_collision_enter
    @on_collision_enter.setter
    def on_collision_enter(self, funs:list[lambda:()]):
        if isinstance(funs, list):
            self._on_collision_enter = funs
        else: self._on_collision_enter = [funs]

    @property
    def colliding_objects(self) -> list:
        return self._colliding_objects
    @colliding_objects.setter
    def colliding_objects(self, objects):
        self._colliding_objects = objects

# GLOBAL METHODS
def update_all(objects:list[PhysicsObject], dt):
    for obj in objects:
        obj.update(dt)
def clear_forces(objects:list[PhysicsObject]):
    for obj in objects:
        obj.clear_force()
def draw_all(objects:list[PhysicsObject], surface):
    for obj in objects:
        obj.draw(surface)
# END GLOBAL METHODS


# Circle
class Circle(PhysicsObject):
    DRAG_C = 0.45
    # for no fill circle, set fill_color to None
    def __init__(self, radius:float=100, fill_color:Color=(255,100,0), outline_color:Color=(255,100,0),  outline_width=1, name="Circle", **kwargs):
        # kwargs = "keyword arguments"
        # ** means pack the unused keyword arguments into a dictionary
        # i.e. mass: 1, pos: (0,0), vel: (0,0)
        self.radius = radius
        self.size = (2*radius, 2*radius)
        super().__init__(name, **kwargs)
        self.fill_color = fill_color if callable(fill_color) else lambda: fill_color
        self.outline_color = outline_color if callable(outline_color) else lambda: outline_color
        self.outline_width = outline_width if callable(outline_width) else lambda: outline_width
        self.contact_type = "Circle"

    def __del__(self):
        super().__del__()

    def update_momi(self):
        self.momi = (self.mass*self.radius**2)/2

    def set_bounds(self):
        self.bounds = Bounds(self.pos-(self.radius,self.radius), (self.radius*2,self.radius*2))

    def update_bounds(self):
        self.bounds.pos = self.pos - (self.radius, self.radius)
        d = 2*self.radius
        self.bounds.size = (d,d)

    def update_shape(self, dt=0):
        return super().update_shape(dt)


    @property # radius
    def radius(self) -> float:
        return self._radius() if callable(self._radius) else self._radius
    @radius.setter
    def radius(self, radius:float):
        self._radius = radius
    
    @property # color
    def fill_color(self) -> Color:
        return self._fill_color() if callable(self._fill_color) else self._fill_color
    @fill_color.setter
    def fill_color(self, color:Color):
        if color is None: color = Color(0,0,0,0)
        self._fill_color = color

    @property # color
    def outline_color(self) -> Color:
        return self._outline_color()
    @outline_color.setter
    def outline_color(self, color:Color):
        c = colors.make_color(color)
        self._outline_color = c if callable(c) else lambda: c
    
    @property # width
    def outline_width(self):
        return self._outline_width()
    @outline_width.setter
    def outline_width(self, width):
        self._outline_width = width if callable(width) else lambda: width
        
    def draw(self, window:Surface):
        radius = max(self.radius, 1)
        s:Surface = Surface((abs(radius*2),abs(radius*2)), pygame.SRCALPHA)

        # pygame.gfxdraw.filled_circle(s, int(self.radius), int(self.radius), int(self.radius), self.fill_color)
        # pygame.gfxdraw.circle(s, int(self.radius), int(self.radius), int(self.radius), self.outline_color)
        pygame.draw.circle(s, self.fill_color, (radius, radius), radius, 0)
        if self.outline_width:
            pygame.draw.circle(s, self.outline_color, (radius, radius), radius, self.outline_width)
        # draw line to see avel
        # pygame.draw.line(s, colors.contrast_color(self.fill_color), (self.radius, self.radius), Vector2(self.radius, self.radius)+Vector2(0,self.radius).rotate(self.angle), 10)
        window.blit(s, self.pos-Vector2(self.radius,self.radius))


# Wall
class Wall(PhysicsObject):
    def __init__(self, point1=(0,0), point2=(0,0), color=colors.white, width=1, name="Wall", normal_length=0):
        self.set_points(point1, point2)  # this also sets self.pos and self.normal
        super().__init__(pos=self.pos, name=name, mass=math.inf)
        self.color = color
        self.width = width
        self.normal_length = normal_length
        self.contact_type = "Wall"

    @property # color
    def color(self):
        return self._color() if callable(self._color) else self._color
    @color.setter
    def color(self, color:Color):
        self._color = color

    def draw(self, window):
        pygame.draw.line(window, self.color, self.point1, self.point2, self.width)
        if self.normal_length > 0:
            pygame.draw.line(window, self.color, self.pos, self.pos + self.normal_length*self.normal, 10) # normal

    def update(self, dt=0):
        super().update(dt)

    def set_points(self, point1=None, point2=None):
        if point1 is not None:
            self.point1 = Vector2(point1)
        if point2 is not None:
            self.point2 = Vector2(point2)
        self.pos = (self.point1 + self.point2)/2
        self.update_normal()

    def update_normal(self):
        self.normal = (self.point2 - self.point1).normalize().rotate(90)

    def set_bounds(self):
        x = min(self.point1.x, self.point2.x)
        y = min(self.point1.y, self.point2.y)
        w = abs(self.point1.x - self.point2.x)
        h = abs(self.point1.y - self.point2.y)
        self.bounds = Bounds((x,y),(w,h))

    def update_bounds(self):
        self.bounds.pos = (min(self.point1.x, self.point2.x), min(self.point1.y, self.point2.y))
        self.bounds.size = (abs(self.point1.x - self.point2.x), abs(self.point1.y - self.point2.y))   
    
    def update_shape(self, dt=0):
        self.point1 += self.vel * dt
        self.point2 += self.vel * dt
        return super().update_shape()

# Polygon
class Polygon(PhysicsObject):
    def __init__(self, local_points=[], fill_color=(255,0,0), outline_color:Color=(255,100,0), outline_width = 1, normals_length=0, **kwargs):
        super().__init__(**kwargs) 
        self.local_points:list[Vector2] = local_points
        self.fill_color = fill_color
        self.outline_color = outline_color
        self.outline_width = outline_width
        self.normals_length = normals_length
        self.contact_type = "Polygon"
        # self.update_polygon() # initialize self.points, NOTE: happens in self.local_points setter

    # PROPERTIES
    @property
    def local_points(self) -> list[Vector2]:
        return self._local_points
    @local_points.setter
    def local_points(self, points):
        self._local_points = [Vector2(p) for p in points]
        self.local_normals:list[Vector2] = []
        for i in range(len(self.local_points)):
            self.local_normals.append((self.local_points[i] - self.local_points[i-1]).normalize().rotate(90))
        self.update_polygon()


    # METHODS
    def update(self, dt=0):
        super().update(dt)

    def draw(self, window):
        if colors.make_color(self.fill_color).a !=0:
            pygame.draw.polygon(window, self.fill_color, self.points, 0)
        if self.outline_width != 0:
            pygame.draw.polygon(window, self.outline_color, self.points, self.outline_width)
        if self.normals_length != 0:
            # for normal, point in zip(self.normals, self.points):
            #     pygame.draw.line(window, self.color, point, point+normal*self.normals_length)
            for i, normal in enumerate(self.normals):
                point = (self.points[i] + self.points[i-1]) /2
                pygame.draw.line(window, self.outline_color, point, point+normal*self.normals_length)

    def check_convex(self):
        if len(self.local_points) > 2:
            n = len(self.local_points)
            convex = True
            for i in range(n):
                # makes a list containing the directionality of [i]->[j] relative to the normal of line [i],[i-1]
                d = [(self.local_points[j%n] - self.local_points[i]).dot(self.local_normals[i])
                     for j in range(i+1, i+n-1)]
                # if each value is <= 0, then there are no points in the normal direction, and angle at point [i] is convex
                if max(d) <= 0:
                    pass
                # otherwise, if each value is >= 0, then the normal direction needs to be flipped because it's pointing the wrong way
                elif min(d) >= 0:
                    self.local_normals[i] *= -1
                # otherwise, if there is a point in the normal direction and a point that is not, the shape must be concave somewhere
                else: convex = False
            if not convex and self.generate_contacts:
                print("WARNING! Non-convex polygon defined. Collisions will be incorrect.")
            return convex
        else:
            print("WARNING! A polygon must have at least 3 points!")
            return False

    def update_polygon(self, check_convex=True):
        if check_convex: self.check_convex()
        # the points of polygon, rotated to object angle
        self.update_points()
        # the current normal, rotated to object angle
        self.update_normals()
        # set the bounding radius (dist from center to furthest point from center)
        self.bounding_radius = 0
        for pt in self.local_points:
            r = pt.magnitude()
            if r > self.bounding_radius: self.bounding_radius = r

    def set_bounds(self):
        return super().set_bounds()

    def update_bounds(self):
        minx = min(point.x for point in self.points)
        maxx = max(point.x for point in self.points)
        miny = min(point.y for point in self.points)
        maxy = max(point.y for point in self.points)
        self.bounds.pos = (minx,miny)
        self.bounds.size = (maxx-minx, maxy-miny)
        
    def update_shape(self, dt=0):
        self.update_polygon(False)
        return super().update_shape(dt)


    # must be called whenever pos, angle, or local_points are changed
    def update_points(self):
        self.points = [p.rotate(self.angle) + self.pos for p in self.local_points]
    
    # must be called whenever angle, or local_points are changed
    def update_normals(self):
        self.normals = [n.rotate(self.angle) for n in self.local_normals]

    def clone(self):
        return copy.deepcopy(self)

    def flip(self, flip_x=True, flip_y=False):
        c = Vector2(-1 if flip_x else 1, -1 if flip_y else 1)
        new_locals = [point*c for point in self.local_points]
        self.local_points = new_locals

    def reflect(self, point1, point2):
        # find the new vertices by getting ortho vector to line
        u:Vector2 = (Vector2(point2) - point1).normalize()
        new_pos = 2*(point1 + ((self.pos-point1)*u)*u) - self.pos
        new_locals = [(2*(point1 + ((point-point1)*u)*u) - point)-new_pos for point in self.points]
        self.pos = new_pos
        self.local_points = new_locals


# Uniform Circle
class UniformCircle(Circle):
    def __init__(self, density=1, radius=100, **kwargs):
        # calculate mass and moment of inertia
        self.density = density
        mass = density * math.pi*radius**2
        momi = 0.5 * mass * radius**2
        if kwargs.get("momi", None) is not None: kwargs.pop("momi")
        super().__init__(radius=radius, mass=mass, momi=momi, **kwargs)

# Uniform Polygon
class UniformPolygon(Polygon):
    def __init__(self, density=None, local_points=[], pos=[0,0], angle=0, shift=True, mass=None, **kwargs):
        if mass is not None and density is not None:
            raise("Cannot specify both mass and density.")
        if density is None:
            density = 1
            if mass is None:
                mass = 1 # if nothing specified, default to mass = 1
        self.density = density
        # print(f"lcl pts: {local_points}")

        # Calculate mass, moment of inertia, and center of mass
        # by looping over all "triangles" of the polygon
        total_mass:float = 0
        total_momi:float = 0
        com_numerator:Vector2 = Vector2(0)
        for i in range(len(local_points)):
            s0 = Vector2(local_points[i])
            s1 = Vector2(local_points[i-1])
            # triangle mass
            tri_area = 0.5 * s0.cross(s1)
            tri_mass =  density * tri_area
            # triangle moment of inertia
            tri_momi = tri_mass/6 * (s0*s0 + s1*s1 + s0*s1)
            # triangle center of mass
            tri_com = (s0 + s1) / 3
            # add to total mass
            total_mass += tri_mass
            # add to total moment of inertia
            total_momi += tri_momi
            # add to center of mass numerator
            com_numerator += tri_mass * tri_com
            
    
        # calculate total center of mass by dividing numerator by denominator (total mass)
        com = com_numerator/total_mass
        # if mass is specified, then scale mass and momi
        if mass is not None:
            total_momi *= mass/total_mass
            total_mass = mass

        # Usually we shift local_points origin to center of mass
        if shift:
            # 1. Shift all local_points by com
            local_points = [p - com for p in local_points]
            # 2. Shift pos
            pos += com.rotate(angle)
            # 3. Use parallel axis theorem to correct the moment of inertia
            #   point_momi = com_momi + mass*(com_to_point.magnitude()**2)
            total_momi -= total_mass * com.magnitude_squared()
                        

        # Then call super().__init__() with those correct values
        super().__init__(mass=abs(total_mass), momi=abs(total_momi), local_points=local_points, pos=pos, angle=angle, **kwargs) 

# # Test UniformPolygon
# shape = UniformPolygon(density=0.01, local_points=[[0,0],[20,0],[20,10],[0,10]])
# print(f"Check mass: {shape.mass} = {0.01*10*20}")  # check mass
# print(f"Check momi: {shape.momi} = {shape.mass/12*(10**2+20**2)}")  # check moment of inertia
# print(shape.local_points) # check if rectangle is centered (checks center of mass)
# print([[-10,-5],[10,-5],[10,5],[-10,5]])  




# -- OTHER -- #
# Particle: circle with a size and color that can change over its lifespan.
class Particle(Circle):
    particles = deque()
    MAX_PARTICLES = 1000
    ltime = pygame.time.get_ticks()
    @classmethod
    def update_all(cls, dt):
        while cls.particles and (len(cls.particles) > cls.MAX_PARTICLES or cls.particles[0].lifetime <= 0):
            cls.particles.popleft()
        for particle in cls.particles:
            particle.update(dt)
        # print(f"num particles: {len(Particle.particles)}")
        # print(f"dt: {pygame.time.get_ticks() - Particle.ltime}")
        cls.ltime = pygame.time.get_ticks()
    @classmethod
    def draw_all(cls, window):
        for p in cls.particles:
            # print("Drawing particle")
            p.draw(window)
    
    def __init__(self, size=10.0, color=(100,120,255,100), lifetime=1000, to_size=1, to_color=(200,100,50,0), name="Particle", mass=1, pos=(0, 0), vel=(0, 0), delta_size=-0.1):
        super().__init__(radius=size/2, name=name, mass=mass, pos=pos, vel=vel)
        self._size = size# if callable(size) else lambda: size
        # self.radius = size/2
        self.color = Color(color)
        self.lifespan = float(abs(lifetime))
        self.lifetime = float(abs(lifetime))
        self._delta_size = delta_size if callable(delta_size) else lambda: (float(to_size - size)/self.lifespan)
        self.curr_delta = self.delta_size
        self.delta_color = ([(Color(to_color)[i] - self.color[i])/self.lifespan for i in range(4)])
        self.lists = []
        Particle.particles.append(self)
        
    def __del__(self):
        super().__del__()
    def draw(self, window):
        s = pygame.Surface((abs(self.size),abs(self.size)), pygame.SRCALPHA)
        color = ([abs(c)%256 for c in self.color])
        pygame.draw.circle(s, color, (self.size/2.0,self.size/2.0), self.size/2.0)
        window.blit(s, self.pos-Vector2(self.size/2.0,self.size/2.0))
    
    @property
    def size(self):
        return float(self._size)
    @size.setter
    def size(self, size):
        self._size = size #if callable(size) else size
    
    @property
    def delta_size(self):
        return self._delta_size()
    @delta_size.setter
    def delta_size(self, delta):
        self._delta_size = delta if callable(delta) else lambda: delta
    
    def update(self, dt):
        # print(f"force: {self.force}")
        self.lifetime -= 1000*dt
        if (self.lifetime <= 0):
            for l in self.lists:
                # print(f"list: {l}")
                l.remove(self)
            # PhysicsObject.phys_objs.remove(self)
            del self
            return
        super().update(dt)
        self.size += self.delta_size * 1000*dt
        self.color = ([self.color[i] + self.delta_color[i] * 1000*dt for i in range(4)])

