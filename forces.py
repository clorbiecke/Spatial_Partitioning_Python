import pygame
from pygame import Surface
from pygame.math import Vector2, Vector3
from physics_objects import  PhysicsObject, Circle
from itertools import chain
import math
import contact

# SingleForce acts only on one object in system. Ex. gravity, air resistance
class SingleForce:
    def __init__(self, objects:list[PhysicsObject]):
        self.objects = objects

    def apply(self):
        # Apply the force to all objects in the objects list
        for obj in self.objects:
            if obj is not None:
                self.force(obj)
    
    @property
    def objects(self) -> list[PhysicsObject] | chain[PhysicsObject]:
        return self._objects()
    @objects.setter
    def objects(self, objects:list[PhysicsObject] | chain[PhysicsObject]):
        self._objects = objects if callable(objects) else lambda: objects

# PairForce is a force that acts between all pairs of objects. Ex. gravity between planets, repulsion, magnetism
class PairForce:
    def __init__(self, objects:list[PhysicsObject]):
        self.objects = objects

    def apply(self):
        # Loop over *all* pairs of objects and apply the calculated force
        # to each object, respecting Newton's 3rd Law.  
        # Use either two nested for loops (taking care to do each pair once)
        for i in range(len(self.objects)-1):
            for j in range(i+1, len(self.objects)):
                self.force(self.objects[i], self.objects[j])
                
        # or use the itertools library (specifically, the function combinations).       

# BondForce  is a force that acts only between specific pairs of objects. Ex. moleculer bonds, spring forces
class BondForce:
    def __init__(self, pairs:list[tuple[PhysicsObject,PhysicsObject]]):
        # pairs has the format 
        # [[obj1, obj2], [obj3, obj4], ... ] 
        # ^ each pair representing a bond between two objects
        if any(len(p) != 2 for p in pairs):
            raise("pairs needs to be a list of pairs")
        self.pairs = pairs

    def apply(self):
        # Loop over all *pairs* from the pairs list.  
        # Apply the force to each member of the pair respecting Newton's 3rd Law.
        for a,b in self.pairs:
            self.force(a,b)
        pass
        
# Add Gravity, SpringForce, SpringRepulsion, AirDrag
class Gravity(SingleForce):
    def __init__(self, gravity:tuple=(0,0), objects:list[PhysicsObject]=[], **kwargs):
        super().__init__(objects=objects, **kwargs)
        self.gravity=Vector2(gravity)

    def force(self, obj: PhysicsObject):
        f = obj.mass * self.gravity
        obj.add_force(f)
        # return f
        
# Force is stiffness * overlap
class Container(SingleForce):
    def __init__(self, pos, width, height, stiffness=100, objects:list[PhysicsObject]=[], **kwargs):
        super().__init__(objects=objects, **kwargs)
        self.pos = Vector2(pos)
        self.width = width
        self.height = height
        self.stiffness = stiffness

    def force(self, obj: PhysicsObject):
        left = self.pos.x - (obj.pos.x - obj.radius) 
        right = self.pos.x + self.width - (obj.pos.x + obj.radius)
        top = self.pos.y - (obj.pos.y - obj.radius)
        bot = self.pos.y + self.height - (obj.pos.y + obj.radius)
        f = Vector2(0)
        if left > 0:
            f += Vector2(left * self.stiffness, 0)
        elif right < 0:
            f += Vector2(right * self.stiffness, 0)

        if top > 0:
            f += Vector2(0, top * self.stiffness)
        elif bot < 0:
            f += Vector2(0, bot * self.stiffness)
        obj.add_force(f)
        # return f
    
    def draw(self, window, color=(255,255,255), width=1):
        pygame.draw.rect(window, color, (self.pos.x, self.pos.y, self.width, self.height), width)
        pass
    
# Potential energy has this form: u(r) = e[(d/r)^2n - (d/r)^n]: (e)psilon is strength, r is distance, d is diameter, n is power
# To get force, we take the gradient(derivative)
class LennardJones(PairForce):
    def __init__(self, strength=1, force_0_dist=1, power=3, objects:list[PhysicsObject]=[], **kwargs):
        super().__init__(objects=objects, **kwargs)
        self.epsilon = strength
        self.sigma = force_0_dist
        self.n = power

    def force(self, a: PhysicsObject, b: PhysicsObject):
        r = a.pos - b.pos
        rmag = r.magnitude() + 10
        d = (self.sigma/rmag)**self.n
        f = self.epsilon*self.n*r.normalize()/rmag*d*(2*d - 1)
        a.add_force(f)
        b.add_force(-f)
        return f
    
    
class HarmonicBonds(BondForce):
    def __init__(self,length=0, stiffness=0, pairs:list[tuple[PhysicsObject,PhysicsObject]]=[], **kwargs):
        super().__init__(pairs=pairs, **kwargs)
        self.length = length
        self.stiffness = stiffness

    def force(self, a: PhysicsObject, b: PhysicsObject):
        r = a.pos - b.pos
        f = -self.stiffness * (r.magnitude() - self.length) * r.normalize()
        a.add_force(f)
        b.add_force(-f)
        # return (f)
    
    def draw(self, window, color=(255,255,255), width=1):
        for a, b in self.pairs:
            pygame.draw.line(window, color, a.pos, b.pos, width)

    
class SpringForce(BondForce):
    def __init__(self,length=0, stiffness=0, damping=0, pairs:list[tuple[PhysicsObject,PhysicsObject]]=[], **kwargs):
        super().__init__(pairs=pairs, **kwargs)
        self.length = length
        self.stiffness = stiffness
        self.damping = damping

    def force(self, a: PhysicsObject, b: PhysicsObject):
        r = a.pos - b.pos
        r_normal = r.normalize() if r.length() != 0 else Vector2(0)
        v = (a.vel - b.vel)*(1.0)
        # f = -self.stiffness * (r.magnitude() - self.length) * r.normalize()
        f = (-self.stiffness * (r.magnitude() - self.length) - (self.damping * v *r_normal)) * r_normal
        a.add_force(f)
        b.add_force(-f)
        # return (f)
    
    def draw(self, window, color=(255,255,255), width=1):
        for a, b in self.pairs:
            o = ((a.pos-b.pos).normalize() * a.radius) if (a.pos-b.pos) != Vector2(0) else Vector2(0)
            x = o.y
            y = o.x
            o = Vector2(x,-y)
            p = Vector2(-x,y)
            pygame.draw.line(window, a.d_color, a.pos + o, b.pos - o, width)
            pygame.draw.line(window, b.d_color, a.pos + p, b.pos - p, width)

class SpringRepulsion(PairForce):
    def __init__(self, strength=1, objects:list[PhysicsObject]=[], **kwargs):
        super().__init__(objects=objects, **kwargs)
        self.strength = strength

    def force(self, a: PhysicsObject, b: PhysicsObject):
        r = a.pos - b.pos
        rmag = r.magnitude()
        R = a.radius + b.radius
        if (rmag > R): return
        f = self.strength * (R-rmag) * r.normalize()
        a.add_force(f)
        b.add_force(-f)
        return f

class AirDrag(SingleForce):
    MAX_WIND_SPD = 6000 #pixels/sec
    def __init__(self, air_density: float=0.0001, wind_vel: tuple=(0,0), objects:list[PhysicsObject]=[], **kwargs):
        super().__init__(objects=objects, **kwargs)
        self.air_density = air_density
        self.wind_vel = Vector2(wind_vel)
        self.inc = 1.0

    def force(self, obj: PhysicsObject):
        v = (obj.vel - self.wind_vel)*(1)
        if v == Vector2(0,0):
            # print(f"{obj.name} vel: {v}")
            return
        f = -0.5 * obj.DRAG_C * self.air_density * (obj.area) * (v*v)
        # -0.5 * 0.5 * 1.224 * math.pi*(obj.radius/270)**2 *
        obj.add_force(f * v.normalize())
        # print(f"{obj.name} vel: {v}, drag force: {f * v.normalize()}")
        # return f
    
    def inc_wind_velx(self):
        self.wind_vel.x += self.inc
    def dec_wind_velx(self):
        self.wind_vel.x -= self.inc
    def set_wind_velx(self, x:float=0):
        self.wind_vel.x = x

    def draw(self, window, bar_pos: tuple=(0,0), bar_max_size: tuple=(100,20), border_color=(0,0,0), fill_color=(0,0,0)):
        # draw bar at top of screen to show wind speed in m/s
        ## draw border
        pygame.draw.rect(window, border_color, (bar_pos, bar_max_size))
        # draw as particles or arrows or something moving at wind_vel on screen
        pass

# Friction force is: F = -v.norm * min(u*m*g, m*v/dt + F*v.norm). u: friction coefficient, m: obj mass, g: force of gravity
# NOTE: must be applied after other forces (uses force on object in calculation)
class Friction(SingleForce):
    def __init__(self, friction:float=0.3, gravity=320*9.8, dt=(1/60000), objects:list[PhysicsObject]=[], **kwargs):
        super().__init__(objects=objects, **kwargs)
        self.friction = friction
        self._gravity = gravity if callable(gravity) else lambda: gravity
        self._dt = dt if callable(dt) else lambda: dt

    def force(self, obj:PhysicsObject):
        if obj.vel == Vector2(0): return 
        u=self.friction
        m=obj.mass
        # g=(self.gravity) # = sqrt(gravity^2 - obj.force.magnitude()^2)
        g = math.sqrt(self.gravity**2 - obj.force.magnitude_squared())
        v=obj.vel
        F=obj.force
        dt=self.dt
        max_friction = (v*m)/dt + F
        f_direction = -v.normalize()#-max_friction.normalize()

        #       usual=> |<--->| |<-------------------->| <=reduced to stop at 0
        F_friction = min(u*m*g, max_friction.magnitude())
        obj.add_force(F_friction*f_direction)
        print(f"Added friction: {F_friction*f_direction} to obj: {obj.name}")
        # print(f"Max friction: {max_friction} calc friction: {u*m*g}")
        # print(f"v after friction: {v + ((F+(F_friction*f_direction))/m)*dt}")
        # to prevent friction force from creating negative velocity...
        # vel after other forces are applied:
        #     vel1 = vel + (force/mass)*dt
        # vel after friction is applied:
        #     vel2 = vel + (force - vel.normal*(f*m*g))*dt
        # max force of friction (f*m*g) such that {vel+((force-(f*m*g))/mass)*dt = 0}
        #     (f*m*g) = -((vel*mass)/dt + force)

    @property
    def dt(self):
        return self._dt()
    @property
    def gravity(self):
        return self._gravity()




#-- OTHER --#
class Bounce:
    def __init__(self, objects:list[PhysicsObject]):
        self.objects = objects

    def apply(self, dt):
        # Loop over *all* pairs of objects and apply the calculated force
        # to each object, respecting Newton's 3rd Law.  
        # Use either two nested for loops (taking care to do each pair once)
        for i in range(len(self.objects)-1):
            for j in range(i+1, len(self.objects)):
                if self.overlapping(self.objects[i], self.objects[j]):
                    print("Overlap detected!")
                    self.resolve(self.objects[i], self.objects[j], dt)
    
    def overlapping(self, o1, o2):
        pass
    def resolve(self, o1, o2, dt):
        pass
class Circle_Collision(Bounce):
    def __init__(self, objects:list[PhysicsObject]):
        self.objects = objects
        # self.object_class = object_class
    
    def overlapping(self, c1:Circle, c2:Circle):
        return (c1.pos.distance_to(c2.pos) < c1.radius + c2.radius)

    def resolve(self, c1:Circle, c2:Circle, dt):
        # Go back in time (record dt subtracted)
        v = c2.vel - c1.vel
        r = c2.pos - c1.pos
        R = c2.radius + c1.radius
        a = (v*v)
        b = (2*r*v)
        c = r*r - R**2
        ## determine how much time to go back to undo overlap (aka, make distance c1->c2 equal R)
        t = (-b + math.sqrt(b**2 - 4*a*c))/(2*a)
        if t > 0:
            t = (-b - math.sqrt(b*b - 4*a*c))/(2*a)
        
        c1.update(t)
        c2.update(t)

        # Apply impulse and update physics with subtracted dt
        n1 = (c2.pos - c1.pos).normalize() # direction of impulse for c1
        n2 = (c1.pos - c2.pos).normalize() # direction of impulse for c2

        ## normal component of velocity: determines how much of momentum is transfered
        nv1 = (c1.vel * n1)/(n1 * n1) * n1
        nv2 = (c2.vel * n2)/(n2 * n2) * n2
        
        ## component of velocity orthogonal to normal component: determins how much of current velocity is retained
        ov1 = c1.vel - nv1
        ov2 = c2.vel - nv2

        M = c1.mass + c2.mass

        c1.vel = ov1 + ((c1.mass - c2.mass)*nv1 + 2*c2.mass*nv2)/M
        c2.vel = ov2 + ((c2.mass - c1.mass)*nv2 + 2*c1.mass*nv1)/M

        c1.update(-t)
        c2.update(-t)

class Slope_Gravity(SingleForce):
    def __init__(self, surface:Surface, gravity:float, objects:list[PhysicsObject]=[], **kwargs):
        super().__init__(objects=objects, **kwargs)
        self.surface:Surface = surface
        self.gravity:float = gravity
        print(f"Gravity: {gravity}")

    def force(self, obj: PhysicsObject):
        # TODO: Alter algorithm to find highest point in radius: get vector from there to current pos
        f = self.gravity
        pos_col = self.surface.get_at([int(obj.pos.x),int(obj.pos.y)])
        xneg_col = self.surface.get_at([int(obj.pos.x)-1,int(obj.pos.y)])
        xpos_col = self.surface.get_at([int(obj.pos.x)+1,int(obj.pos.y)])
        yneg_col = self.surface.get_at([int(obj.pos.x),int(obj.pos.y)-1])
        ypos_col = self.surface.get_at([int(obj.pos.x),int(obj.pos.y)+1])
        z_x = -((xpos_col.a * int(xpos_col.r/255))-(xneg_col.a * int(xneg_col.r/255)))
        z_y = -((ypos_col.a * int(ypos_col.r/255))-(yneg_col.a * int(yneg_col.r/255)))
        Fx = f * z_x/Vector2(2,z_x).magnitude()
        Fy = f * z_y/Vector2(2,z_y).magnitude()
        slope = Vector3(2,2,(z_x+z_y)/2).normalize()
        F = (f * slope.z / Vector2(slope.x, slope.y).magnitude()) *  Vector2(slope.x, slope.y)
        # obj.add_force(Vector2(Fx, Fy))
        obj.add_force(F)
        print(f"Added grav: {F}")
        # print(f"Added grav: {Vector2(Fx, Fy)}")