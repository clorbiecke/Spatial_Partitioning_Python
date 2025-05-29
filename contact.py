from pygame.math import Vector2
from physics_objects import Circle, Wall, PhysicsObject, Polygon
import math

# Returns a new contact object of the correct subtype
# This function has been done for you.
def generate(a, b, addl_mod:str="",**kwargs) -> 'Contact':
    # Check if a's type comes later than b's alphabetically.
    # We will label our collision types in alphabetical order, 
    # so the lower one needs to go first.
    if not a.generate_contacts or not b.generate_contacts: return None
    if b.contact_type < a.contact_type:
        a, b = b, a
    # This calls the class of the appropriate name based on the two contact types.
    return globals()[f"{a.contact_type}_{b.contact_type}{f"_{addl_mod}" if addl_mod != "" else ""}"](a, b, **kwargs)
    
# Generic contact class, to be overridden by specific scenarios
class Contact():
    def __init__(self, a:PhysicsObject, b:PhysicsObject, resolve=False, **kwargs):
        self.a:PhysicsObject = a
        self.b:PhysicsObject = b
        self.kwargs = kwargs
        self.on_impact = []
        self.update()
        self.bool = self.overlap > 0
        if self.bool: 
            self.collision_enter()
            # print(f"{a.name}({a.__class__}) collided with {b.name}({b.__class__})")
        if resolve:
            # print(f"self.bool = {bool(self.resolve(update=False))}")
            self.bool = bool(self.resolve(update=False))

    @property
    def on_impact(self):
        return self._on_impact
    @on_impact.setter
    def on_impact(self, func:list[lambda:None]):
        self._on_impact = func
    def impact(self, **kwargs):
        for f in self.on_impact:
            f(**kwargs)
    
    def __bool__(self):
        return self.bool

    def update(self):  # virtual function
        self.overlap = 0
        self.normal = Vector2(0, 0)

    def resolve(self, restitution=None, rebound=None, friction=None, update=True):
        if update:
            self.update()
        restitution = restitution if restitution is not None else self.kwargs.get("restitution", 0)
        rebound = rebound if rebound is not None else self.kwargs.get("rebound", 0) 
        friction = friction if friction is not None else self.kwargs.get("friction", 0) 
        # ^ priority for restitution is: 1 argument in resolve, 2 argument in generate, 3 default value = 1 (elastic)
        if self.overlap > 0:
            a:PhysicsObject = self.a
            b:PhysicsObject = self.b
            # resolve overlap
            self.undo_overlap()
            # resolve velocity: make sure movement is toward eachother
            return self.resolve_velocity(restitution, rebound, friction)
        return False
    
    # helper for resolve()
    def undo_overlap(self):
        # print("Contact.undo_overlap()")
        a:PhysicsObject = self.a
        b:PhysicsObject = self.b
        # m = 1/((1/a.mass)+(1/b.mass))
        # o = self.normal * self.overlap
        # a.delta_pos(o * m/a.mass)
        # b.delta_pos(-o * m/b.mass)
        # RESOLVE OVERLAP
        m = 1/(1/self.a.mass + 1/self.b.mass) # reduced mass
        self.a.set(pos=self.a.pos + m/self.a.mass*self.overlap*self.normal)
        self.b.set(pos=self.b.pos - m/self.b.mass*self.overlap*self.normal)

    # helper for resolve()
    def resolve_velocity(self, restitution, rebound, friction):
        a:PhysicsObject = self.a
        b:PhysicsObject = self.b

        n = self.normal # direction to 'a' from 'b'
        # m = 1/((1/a.mass)+(1/b.mass)) # Reduced Mass (just helpful for calculations)
        # Reduced mass for collisions with rotation
        Sa:Vector2 = self.point()-self.a.pos
        Sb:Vector2 = self.point()-self.b.pos
        m = 1/((1/a.mass) + (1/b.mass) + (Sa.cross(n)**2 /a.momi) + (Sb.cross(n)**2 /b.momi))
        # calculate relative velocity at contact point
        vpt_a = math.radians(a.avel)*Vector2(self.point()-self.a.pos).rotate(90)
        va = a.vel + vpt_a
        vpt_b =  math.radians(b.avel)*Vector2(self.point()-self.b.pos).rotate(90)
        vb = b.vel + vpt_b
        v = va - vb # vel of 'a' relative to 'b'

        # v = a.vel - b.vel # vel of 'a' relative to 'b'
        if v * n < 0:
            # print("resolving velocity...")
            # Calculate impulse to apply to objects for bounce
            Jn = -(1 + restitution)*m*v*n + m*rebound
            impulse = Jn*n
            #NOTE: commented out because friction does not work with angular impulse calculation
            if friction != 0:
                # Calculate impulse to apply for friction
                # (1) Calculate tangential normal, make sure v dot t is <= 0
                t = self.normal.rotate(90)
                if v * t > 0: t = -t
                # (2) Calculate impulse needed to make vf dot t == 0
                Jf = -friction*v*t
                # (3) Check if static or kinetic is neccessary
                if Jf <= friction * Jn: # use static
                    # If there is static friction, correct for creep
                    correction = ((v*t)/(v*n))*self.overlap
                    self.a.pos += m/self.a.mass * correction * t
                    self.a.pos += m/self.b.mass * correction * t
                else: 
                    Jf = friction * Jn

                # (4) Calculate the total impulse
                impulse += Jf*t

            # apply impulses
            a.add_impulse(impulse, self.point())
            # print(f"Added impulse {J*n} to circle.")
            b.add_impulse(-impulse, self.point())
            # return True if impulses need to be applied, return False otherwise
            # print(f"impulse: {impulse}")
            return impulse
        else: return False
    
    # call collision_enter() on both objects
    def collision_enter(self):
        self.a.collision_enter(other=self.b, contact=self)
        self.b.collision_enter(other=self.a, contact=self)

    # virtual function to return point of contact
    def point(): return
#-- CIRCLE-CIRCLE --#
# Contact class for two circles
class Circle_Circle(Contact):
    def update(self):  # compute the appropriate values
        a:Circle = self.a
        b:Circle = self.b
        self.overlap = (a.radius + b.radius) - a.pos.distance_to(b.pos)
        self.normal = (a.pos - b.pos).normalize() if (a.pos - b.pos) != Vector2(0) else Vector2(1,0) # direction to a.pos from contact position (b.pos)

    # def resolve(self, restitution=None, update=True):
    #     return super().resolve(restitution=restitution, update=update)
    
    def point(self):
        return (self.a.pos - self.a.radius*self.normal)

# Complex Circle_Circle collision resolution that involves going back in time to point of contact
class Circle_Circle_Complex(Circle_Circle):
    def update(self):  # compute the appropriate values
        super().update()

    # Overrides Contact.resolve()
    def resolve(self, restitution=None, update=True):
        if update:
            self.update()
        restitution = restitution if restitution is not None else self.kwargs.get("restitution", 0) 
        # ^ priority for restitution is: 1 argument in resolve, 2 argument in generate, 3 default value = 1 (elastic)
        J_applied = False
        if self.overlap > 0:
            a:PhysicsObject = self.a
            b:PhysicsObject = self.b 
            # Undo overlap
            dt = self.undo_overlap() # also updates self.normal
            if dt == 0: return super().resolve(restitution=restitution, update=update)
            # Add impulse
            m = 1/((1/a.mass)+(1/b.mass)) # Reduced Mass: == (a.mass * b.mass)/(a.mass + b.mass)
            ## resolve velocity: make sure movement is toward eachother
            v = a.vel - b.vel # vel of 'a' relative to 'b'
            # if v == Vector2(0): return super().resolve(restitution=restitution, update=update)
            n = self.normal # direction to 'a' from 'b'
            if v * n < 0:
                J = -(1 + restitution)*m*v*n
                # apply impulses
                a.add_impulse(J*n)
                b.add_impulse(-J*n)
                npos = a.pos + (b.pos - a.pos)*(a.radius/(a.radius+b.radius))
                # print(f"J: {J}, npos: {npos}, a.pos: {a.pos}")
                self.impact(impact_normal=self.normal, impact_pos=npos, impact_strength=J)
                J_applied = True
            # Update objects with t removed in undo_overlap
            a.update(-dt)
            b.update(-dt)
        return J_applied
    
    # Go back in time to undo overlap (returns dt)
    def undo_overlap(self):
        c1:Circle = self.a
        c2:Circle = self.b
        v = c1.vel - c2.vel # vel of c1 relative to c2
        if v == Vector2(0): return 0
        r = c1.pos - c2.pos # pos of c1 relative to c2
        R = c1.radius + c2.radius # distance we want between c1 and c2
        ## coefficients in 0 = (v^2)t^2 + (2*r*v)t + (r*r - R**2)
        a = (v*v) 
        b = (2*r*v)
        c = r*r - R**2
        ## determine how much time to go back to undo overlap (aka, make distance c1->c2 equal R)
        ## calculate 't' using quadratic formula. Use most negative 't'.
        t=1
        try:
            t = (-b + math.sqrt(b**2 - 4*a*c))/(2*a)
        except (ValueError, ZeroDivisionError):
            t=1
            print(f"Error: ({c1.name}, {c2.name}). a: {a}, b: {b}, c: {c}")
            print(f"Trying with '-'...")
        if t > 0:
            try:
                t = (-b - math.sqrt(b**2 - 4*a*c))/(2*a)
            except (ValueError, ZeroDivisionError):
                print(f"Error: ({c1.name}, {c2.name}). a: {a}, b: {b}, c: {c}")
                print(f"\treturning 0 and resolving with less complex algorithm...")
                return 0
                # raise ValueError(f"Error: ({cir.name}, {wall.name}). a: {a}, b: {b}, c: {c}")
        # if t > 0:
        #     t = (-b - math.sqrt(b*b - 4*a*c))/(2*a)
        ## go back in time
        c1.update(t)
        c2.update(t)
        self.normal = (c1.pos - c2.pos).normalize() # update normal
        ## return dt (dt = t)
        return t
#-------------------#


#-- CIRCLE-WALL --#
# Contact class for Circle and an infinitly wide and deep Wall
# Circle is before Wall because it comes before it in the alphabet
class Circle_Wall(Contact):
    def update(self):  # compute the appropriate values
        cir:Circle = self.a
        wall:Wall = self.b
        self.normal = wall.normal
        self.overlap = cir.radius - (cir.pos - wall.point1)*wall.normal
    
    def point(self):
        return (self.a.pos - self.a.radius*self.normal)
        
# Complex Circle_Wall collision resolution that involves going back in time to point of contact
class Circle_Wall_Complex(Circle_Wall):
    def update(self):  # compute the appropriate values
        super().update()

    def resolve(self, restitution=None, friction=None, update=True):
        # FIXME: does not work with friction...
        if update:
            self.update()
        restitution = restitution if restitution is not None else self.kwargs.get("restitution", 0)
        friction = friction if friction is not None else self.kwargs.get("friction", 0) 
        # ^ priority for restitution is: 1 argument in resolve, 2 argument in generate, 3 default value = 1 (elastic)
        J_applied = False
        if self.overlap > 0:
            a:Circle = self.a
            b:Wall = self.b 
            # Undo overlap
            dt = self.undo_overlap() # also updates self.normal
            # undo_overlap() returns 0 if there is an issue (like velocity 0)
            if dt == 0: 
                print("dt was 0, returning super.resolve()")
                return super().resolve(restitution=restitution, friction=friction, update=update)
            else:
                J_applied = self.resolve_velocity(restitution, friction)
                if J_applied:
                    self.impact(impact_normal=self.normal, impact_pos=a.pos-(self.normal*a.radius), impact_strength=J_applied)
                # Update objects with -dt from undo_overlap to bring back to present time
                a.update(-dt)
                b.update(-dt)
        return J_applied
    
    def undo_overlap(self):
        cir:Circle = self.a
        wall:Wall = self.b
        t = 0
        v = cir.vel - wall.vel # vel of cir relative to wall
        if v == Vector2(0): 
            # calculation impossible if relative velocity is 0
            return 0 
        # r = cir.pos - self.npos # pos of cir relative to wall
        # R = cir.radius # distance we want between cir and wall
        # Undo to when circle touches point on line end
        # if self.npos == wall.point1 or self.npos == wall.point2:
        #     ## coefficients in 0 = (v^2)t^2 + (2*r*v)t + (r*r - R**2)
        #     a = (v*v) 
        #     b = (2*r*v)
        #     c = r*r - R**2
        #     ## determine how much time to go back to undo overlap (aka, make distance cir->wall equal R)
        #     ## calculate 't' using quadratic formula. Use most negative 't'.
        #     t=1
        #     try:
        #         t = (-b + math.sqrt(b**2 - 4*a*c))/(2*a)
        #     except (ValueError, ZeroDivisionError):
        #         t=1
        #         print(f"Error: ({cir.name}, {wall.name}). a: {a}, b: {b}, c: {c}")
        #         print(f"Trying with '-'...")
        #     if t > 0:
        #         try:
        #             t = (-b - math.sqrt(b**2 - 4*a*c))/(2*a)
        #         except (ValueError, ZeroDivisionError):
        #             print(f"Error: ({cir.name}, {wall.name}). a: {a}, b: {b}, c: {c}")
        #             print(f"\treturning 0 and resolving with less complex algorithm...")
        #             return 0
        #             # raise ValueError(f"Error: ({cir.name}, {wall.name}). a: {a}, b: {b}, c: {c}")

        #     cir.update(t)
        #     wall.update(t)
        #     ## check if still overlapping
        #     self.update()
        #     if self.overlap <= 0: return t
        #     else:
        #         print(f"ERROR: overlap resolution between {cir.name} and {wall.name} failed. Using generic method...")
        #         cir.update(t)
        #         wall.update(t)
        #         self.update()
        #         return 0
        # # if contact pos is on line...
        # else:
        # Undo to when circle touches point on line
        d = (self.normal * self.overlap)
        t = (d*d)/(v*d)
        cir.update(t)
        wall.update(t)
        self.update()
        ## return dt (dt = t)
        return t
#-----------------#


#-- CIRCLE-POLYGON --#
# Contact class for circle and polygon
class Circle_Polygon(Contact):
    def update(self):
        # print("circle_poly collision")
        cir:Circle = self.a
        poly:Polygon = self.b
        self.overlap = math.inf
        side = 0
        for i in range(len(poly.points)):
            # check for overlap
            overlap = cir.radius - (cir.pos - poly.points[i])*poly.normals[i]
            # if overlapping, determine if its on line or near a point
            if overlap < self.overlap:
                side = i
                self.overlap = overlap
                self.normal = poly.normals[i]
        # in class version:
        if self.overlap > 0 and self.overlap < cir.radius:
            # rather than project onto line, we check if the angle between p2->p1 and p1->cir is less than 90 degrees
            # use dot product. if > 0, then circle is beyond endpoint
            if (poly.points[side] - poly.points[side-1]) * (cir.pos - poly.points[side]) > 0:
                self.normal = (cir.pos - poly.points[side]).normalize()
                self.overlap = cir.radius - poly.points[side].distance_to(cir.pos)
            elif (poly.points[side-1] - poly.points[side]) * (cir.pos - poly.points[side-1]) > 0:
                self.normal = (cir.pos - poly.points[side-1]).normalize()
                self.overlap = cir.radius - poly.points[side-1].distance_to(cir.pos)


    def point(self):
        return (self.a.pos - self.a.radius*self.normal)

# Complex Circle_Polygon collision resolution that involves going back in time to point of contact
class Circle_Polygon_Complex(Circle_Polygon):
    def update(self):  # compute the appropriate values
        super().update()

    def resolve(self, restitution=None, update=True):
        if update:
            self.update()
        restitution = restitution if restitution is not None else self.kwargs.get("restitution", 0) 
        # ^ priority for restitution is: 1 argument in resolve, 2 argument in generate, 3 default value = 1 (elastic)
        J_applied = False
        if self.overlap > 0:
            cir:Circle = self.a
            poly:Polygon = self.b
            # Undo overlap
            dt = self.undo_overlap() # also updates self.normal
            # undo_overlap() returns 0 if there is an issue (like velocity 0)
            if dt == 0: 
                return super().resolve(restitution=restitution, update=update)
            else:
                J_applied = self.resolve_velocity(restitution)
                if J_applied:
                    self.impact(impact_normal=self.normal, impact_pos=cir.pos-(self.normal*cir.radius), impact_strength=J_applied)
                # Update objects with -dt from undo_overlap to bring back to present time
                cir.update(-dt)
                poly.update(-dt)
        return J_applied

    def undo_overlap(self):
        cir:Circle = self.a
        poly:Polygon = self.b
        v = cir.vel - poly.vel
        if v == Vector2(0):
            return 0
        
        t = -math.inf
        n = v.normalize()
        cir1 = Vector2(-n.y, n.x) * cir.radius
        cir2 = Vector2(n.y, -n.x) * cir.radius
        cirs = [cir1, cir2]
        lines = []
        contact_index = None
        # we want the largest dt (least negative)

        # determine which lines on poly the circle should have collided with
        # for i in range(len(poly.points)):
        #     # velocity of cir * normal should be neg
        #     # TODO: fix to make sure to check for collisions with both endpoints of each line
        #     if poly.normals[i] * v < 0:
        #         # check if line overlaps circle vector line
        #         p_p2 = poly.points[i-1]-poly.points[i]
        #         p_cir1 = cir1 - poly.points[i]
        #         p_cir2 = cir2 - poly.points[i]
        #         ## proj pos onto line
        #         p_cir1_proj = (p_cir1 * p_p2)/(p_p2 * p_p2)
        #         p_cir2_proj = (p_cir2 * p_p2)/(p_p2 * p_p2)
        #         ## get line component of pos vector of cir rel to p
        #         proj_comp1 = (p_cir1_proj * p_p2)
        #         proj_comp2 = (p_cir2_proj * p_p2)
        #         ## get ortho component of pos vector of cir rel to p
        #         ortho_comp1 = p_cir1 - (p_cir1_proj * p_p2)
        #         ortho_comp2 = p_cir2 - (p_cir2_proj * p_p2)
        #         ## get pos of cir proj onto line
        #         proj_pos1 = poly.points[i] + (p_cir1_proj * p_p2)
        #         proj_pos2 = poly.points[i] + (p_cir2_proj * p_p2)
        #         ## use proj_pos and vel to calculate time to intersection
        #         ### t = (cir->proj_pos).length / (vel proj onto cir->proj_pos)
        #         c1_pp1 = (proj_pos1 - cir1)
        #         v1 = (v*c1_pp1)/(c1_pp1*c1_pp1)
        #         t1 = 1 / v1
        #         int_cir1 = cir1 + v*t1

        #         c2_pp2 = (proj_pos2 - cir2)
        #         v2 = (v*c2_pp2)/(c2_pp2*c2_pp2)
        #         t2 = 1 / v2
        #         int_cir2 = cir2 + v*t2
        #         ## See if cir will intersect line
        #         num_int = 0
        #         if int_cir1 - poly.points[i] > 0: num_int += 1
        #         if int_cir2 - poly.points[i] > 0: num_int += 1
        #         if num_int > 0:
        #             # if furthest point on cir ortho to line intersects line then use it as contact point
        #             # otherwise, determine dist to collision with (vel comp of (cir->point)).len + sqrt(rad^2 - (ortho comp).len^2)
        #             overlap_dist = cir.radius - (cir.pos - poly.points[i])*poly.normals[i]
        #             dt = (1/((v*(poly.normals[i]*overlap_dist))/((poly.normals[i]*overlap_dist)*(poly.normals[i]*overlap_dist))))
        #             int_pos = cir.pos + v * dt
        #             if (int_pos - poly.points[i]) * p_p2.normalize() > 0:
        #                 # then contact point should be furthest point on cir ortho to line
        #                 if dt > t and dt < 0:
        #                     t = dt  
        #                     contact_point = i
  
        #             else:
        #                 # TODO: determine dist to collision with (vel comp of (cir->point)).len + sqrt(rad^2 - (ortho comp).len^2)
        #                 #       and determine dt with that distance
        #                 vel_comp = (((poly.points[i] - cir.pos)*v)/(v*v))*v
        #                 overlap_dist = vel_comp.length() + math.sqrt(cir.radius**2 - ((poly.points[i] - cir.pos)-vel_comp).length_squared())
        #                 dt = -overlap_dist/v.length
        #                 if dt > t and dt < 0:
        #                     t = dt
        #                     contact_index = i


        # FIXME: IN PROGRESS
        for i in range(len(poly.points)):
            # TODO: fix to make sure to check for collisions with both endpoints of each line
            # velocity of cir * normal should be neg
            if poly.normals[i] * v < 0:
                # IDEA: start by getting furthest point of cir that is ortho to line and going back in time to put it on the line
                #       then, if it is actually on the line (between two points) use that point to find dt
                #       otherwise, find component of (cir.pos->nearest point) ortho to vel
                #       (if comp.len > radius, then there is no collision, you can skip the rest)
                #       then use radius - sqrt(radius^2 - comp.len^2) to get distance between the circle's nearest point 
                #       and the line point along vel, then divide that dist by v.len and add that to dt from prev step
                # proj p->cir.pos onto p->p2 and subtract from p->cir.pos to get ortho comp
                # FIXME: replace cir.pos below with p + (ortho comp + (radius * ortho comp normalized)
                # divide length of the ortho comp by scalar of v proj onto ortho comp
                # multiply that by v and add cir.pos to get pos of cir intersecting line

                p_cir = cir.pos - poly.points[i]
                p_cir_ortho_line = ((p_cir*-poly.normals[i])/(poly.normals[i]*poly.normals[i])) * -poly.normals[i]
                far_point_ortho_line = (p_cir_ortho_line * (1 + (cir.radius/p_cir_ortho_line.magnitude())))
                1/((v * far_point_ortho_line) / (far_point_ortho_line*far_point_ortho_line))
                # IN PROGRESS



        # go back in time
        cir.update(t)
        poly.update(t)
        # alter normal and overlap to use that line
        self.update()

        # return t
        return t
#--------------------#            


#-- WALL-WALL --#
# Empty class for Wall - Wall collisions
# The intersection of two infinite walls is not interesting, so skip them
class Wall_Wall(Contact):
    pass
class Wall_Wall_Complex(Wall_Wall):
    pass
#---------------#


#-- POLYGON-WALL --#
# Empty class for Polygon - Wall collisions
class Polygon_Wall(Contact):
    # calculate/set self.overlap and self.normal
    def update(self):
        poly:Polygon = self.a
        wall:Wall = self.b
        self.overlap = -math.inf
        self.normal = wall.normal
        self.index = None
        for i in range(len(poly.points)):
            r = -(poly.points[i] - wall.point1)*wall.normal
            if r > self.overlap: 
                self.overlap = r
                self.index = i


    # get the point of contact
    def point(self):
        poly:Polygon = self.a
        return poly.points[self.index]

class Polygon_Wall_Complex(Polygon_Wall):
    pass
#---------------#

#-- POLYGON-POLYGON --#
# Empty class for Polygon - Polygon collisions
class Polygon_Polygon(Contact):
    # calculate/set self.overlap and self.normal
    def update(self):
        a:Polygon = self.a
        b:Polygon = self.b
        self.overlap = math.inf
        self.normal = None
        self.index = None
        self.polygon = None

        # for each side in 'b': calculate overlap with poly 'a'
        # if overlap is less than self.overlap, save it

        # for each wall in 'b'
        for j in range(len(b.points)):
            overlap = -math.inf
            index = None
            # find point in 'a' with most overlap
            for i in range(len(a.points)):
                r = -(a.points[i] - b.points[j])*b.normals[j]
                if r > overlap:
                    overlap = r
                    index = i
            # if overlap < self.overlap: save
            if overlap < self.overlap:
                self.overlap = overlap
                self.normal = b.normals[j]
                self.index = index 
                self.polygon = a
            # if overlap <= 0, that means all points in 'a' must be in normal direction, so overlap is impossible
            if overlap <= 0: return False
        # flip roles: treat 'b' as poly and 'a' as walls
        for j in range(len(a.points)):
            overlap = -math.inf
            index = None
            # find point in 'a' with most overlap
            for i in range(len(b.points)):
                r = -(b.points[i] - a.points[j])*a.normals[j]
                if r > overlap:
                    overlap = r
                    index = i
            # if overlap < self.overlap: save
            if overlap < self.overlap:
                self.overlap = overlap
                self.normal = -a.normals[j]
                self.index = index 
                self.polygon = b
        # print(f"Overlap: {self.overlap}")   

            
    # get the point of contact
    def point(self):
        return self.polygon.points[self.index]

class Polygon_Polygon_Complex(Polygon_Polygon):
    pass
#---------------#

class Ray:
    def __init__(self, origin:Vector2, direction:Vector2, magnitude=10000):
        self.origin = origin
        self.direction = direction
        self.magnitude = magnitude
        self.contact_type = "Ray"

class RayCastHit:
    def __init__(self, object, ray, hit_pos, to_hit):
        self.object:PhysicsObject = object
        self.ray:Ray = ray
        self.hit_pos:Vector2 = hit_pos
        self.to_hit:Vector2 = to_hit

class RayCast:
    def __init__(self, ray, other):
        self.hit:RayCastHit = None
        self.ray = ray
        self.other = other
        self.update()

    def update(self):
        pass

def raycast(other, origin, direction, magnitude=10000):
    ray = Ray(origin, direction, magnitude)
    a = other.contact_type
    b = ray.contact_type

    if b < a:
        a, b = b, a
    # This calls the class of the appropriate name based on the two contact types.
    return globals()[f"{a}_{b}"](ray, other)

class Polygon_Ray(RayCast):
    def update(self):
        poly:Polygon = self.other
        ray:Ray = self.ray

        self.dist = math.inf

        rn = ray.direction.rotate(90).normalize()
        
        for i in range(len(poly.points)):
            #project onto ray normal
            p1 = (poly.points[i] - ray.origin) * rn
            p2 = (poly.points[i-1] - ray.origin) * rn
            if p1*p2 >= 0: continue
            # otherwise, the ray overlaps the poly edge somewhere along its line (not neccessarily the ray)
            
            #project ray onto poly
            r1 = (ray.origin - poly.points[i]) * poly.normals[i]
            r2 = (ray.origin+(ray.direction*ray.magnitude) - poly.points[i]) * poly.normals[i]
            if r1*r2 >= 0:continue

            #calculate distance to collision
            dist = ray.magnitude * (r1/(r1-r2))
            if dist < self.dist and dist > 0: 
                self.dist = dist
                edge = i
        
        if self.dist is not math.inf and self.dist > 0:
            to_hit = (ray.direction*self.dist)
            self.hit = RayCastHit(poly, ray, ray.origin+to_hit, to_hit=to_hit)


class Ray_Wall(RayCast):
    def update(self):
        wall:Wall = self.other
        ray:Ray = self.ray

        self.dist = ray.magnitude

        rn = ray.direction.rotate(90).normalize()
        #project onto ray normal
        p1 = (wall.point1 - ray.origin) * rn
        p2 = (wall.point2 - ray.origin) * rn
        if p1*p2 >= 0: return
        # otherwise, the ray overlaps the poly edge somewhere along its line (not neccessarily the ray)
        
        #project ray onto poly
        r1 = (ray.origin - wall.point1) * wall.normal
        r2 = (ray.origin+(ray.direction*ray.magnitude) - wall.point1) * wall.normal
        if r1*r2 >= 0: return

        #calculate distance to collision
        dist = ray.magnitude * (r1/(r1-r2))
        if dist < self.dist and dist > 0: 
            self.dist = dist
        
        if self.dist is not math.inf and self.dist > 0:
            to_hit = (ray.direction*self.dist)
            self.hit = RayCastHit(wall, ray, ray.origin+to_hit, to_hit=to_hit)




#--------Find Collision Position--------#
# c2->c1: sqrt((c2.x - c1.x)^2 + (c2.y - c1.y)^2)
# we want to know what 't' to multiply velocities by to make the distance = c1.radius + c2.radius = R
# c2 relative velocity: c2.vel - c1.vel
# c2 relative position: c2.pos - c1.pos
# c2->c1 new distance: sqrt((c2.x+(c2.vel.x*t) - c1.x)^2 + (c2.y+(c2.vel.x*t) - c1.y)^2) = R
##   R^2 = (c2.x+(c2.vel.x*t) - c1.x)^2 + (c2.y+(c2.vel.x*t) - c1.y)^2
##   r  = c2.pos - c1.pos  # c2 relative pos
##   v  = c2.vel - c1.vel  # c2 relative vel
##   rx = c2.x - c1.x
##   ry = c2.y - c1.y
##   R^2 = (rx + c2.vel.x*t)^2 + (ry + c2.vel.y*t)^2
##   Distribute: ((rx)^2 + (c2.vel.x)^2 * t^2 + 2(rx*c2.vel.x) * t) + ((ry)^2 + (c2.vel.y)^2 * t^2 + 2(ry*c2.vel.y) * t) - R^2
##   Simplify: [(c2.vel.x)^2 + (c2.vel.y)^2]*t^2 + [2(rx*c2.vel.x) + 2(ry*c2.vel.y)]*t + [(rx)^2 + (ry)^2 - R^2]
##   Even Simpler: (v^2)*t^2 + (2*r*v)*t + (r^2 - R^2) = 0
##                  a^        b^              c^
# Solve for 't': t = (-b +- sqrt(b^2 - 4*a*c))/(2*a)
# c2 final position: c2.pos + c2.vel*t

## in code ##
# c1:Circle = self.a
# c2:Circle = self.b
# # Go back in time (record dt subtracted)
# v = c2.vel - c1.vel
# r = c2.pos - c1.pos
# R = c2.radius + c1.radius
# a = (v*v)
# b = (2*r*v)
# c = r*r - R**2
# ## determine how much time to go back to undo overlap (aka, make distance c1->c2 equal R)
# t = (-b + math.sqrt(b**2 - 4*a*c))/(2*a)
# if t > 0:
#     t = (-b - math.sqrt(b*b - 4*a*c))/(2*a)

# c1.update(t)
# c2.update(t)

# # Apply impulse and update physics with subtracted dt
# n1 = (c2.pos - c1.pos).normalize() # direction of impulse for c1
# n2 = (c1.pos - c2.pos).normalize() # direction of impulse for c2

# ## normal component of velocity: determines how much of momentum is transfered
# nv1 = (c1.vel * n1)/(n1 * n1) * n1
# nv2 = (c2.vel * n2)/(n2 * n2) * n2

# ## component of velocity orthogonal to normal component: determins how much of current velocity is retained
# ov1 = c1.vel - nv1
# ov2 = c2.vel - nv2

# M = c1.mass + c2.mass

# c1.vel = ov1 + ((c1.mass - c2.mass)*nv1 + 2*c2.mass*nv2)/M
# c2.vel = ov2 + ((c2.mass - c1.mass)*nv2 + 2*c1.mass*nv1)/M

# c1.update(-t)
# c2.update(-t)
