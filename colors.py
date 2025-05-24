from __future__ import annotations
import pygame

class Color(pygame.Color):
    def __init__(self, *args):
        return super().__init__(*args)
    
    # Makes the color lighter by increasing color value by (255-value)*norm_percent. 
    # Ex. if norm_percent = 0.5 and r = 155, this function will make r = 205.
    def lighten(self, norm_percent:float=0.5):
        if norm_percent > 0 and norm_percent <= 1:
            self.r = int((255-self.r) * norm_percent)
            self.g = int((255-self.g) * norm_percent)
            self.b = int((255-self.b) * norm_percent)
        else: print(f"lighten() requires a norm_percent value in range (0.0 - 1.0). Given value: {norm_percent}")
        return self

    # Makes the color darker by decreasing color value by value*norm_percent. 
    # Ex. if norm_percent = 0.5 and r = 200, this function will make r = 100.
    def darken(self, norm_percent:float=0.5):
        if norm_percent > 0 and norm_percent <= 1:
            self.r = int(self.r*norm_percent)
            self.g = int(self.g*norm_percent)
            self.b = int(self.b*norm_percent)
        else: print(f"darken() requires a norm_percent value in range (0.0 - 1.0). Given value: {norm_percent}")
        return self
    
    def contrast_lighten_darken(self, norm_percent:float=0.5):
        lumin = self.luminance()
        if lumin > 0.5:
            return self.darken(norm_percent)
        else: return self.lighten(norm_percent)

    # Returns the luminance of this color, a float in range (0.0 - 1.0).
    def luminance(self) -> float:
        r,g,b,_ = self.normalize()
        lumin = 0.2126*r + 0.7152*g + 0.0722*b
        return lumin


clear = (0,0,0,0)
white = (255,255,255)
grey = (160,160,160)
black = (0,0,0)
brown = (165,42,42)

red = (255,0,0)
orange = (255,165,0)
yellow = (255,255,0)
green = (0,255,0)
blue = (0,0,255)
purple = (128,0,128)

cyan = (0,255,255)
majenta = (255,0,255)
pink = (255,192,203)

grass = (0, 190, 0)
cream = (255,253,208)
midnight = (25,25,112)


#-- COLOR EFFECTS --#
def lighten(color:tuple, norm_percent:float=0.5) -> Color | False:
    c = validate_color(color)
    if c is False: return False
    
    return c.lighten(norm_percent)

def darken(color:tuple, norm_percent:float=0.5) -> Color | False:
    c = validate_color(color)
    if c is False: return False
    
    return c.darken(norm_percent)

def contrast_lighten_darken(color:tuple, norm_percent:float=0.5) -> Color | False:
    c = validate_color(color)
    if c is False: return False
    
    return c.contrast_lighten_darken(norm_percent)

def dark(color:Color) -> Color:
    c:Color = Color(make_color(color))
    for i in range(3):
        c[i] = int(c[i]/2)%256

    return c

def light(color:Color) -> Color:
    c:Color = Color(make_color(color))
    for i in range(3):
        c[i] += int((255 - c[i])/2)

    return c

def alpha(color:Color, alpha:int) -> Color:
    c:Color = Color(make_color(color))
    a = int(alpha)
    c.a = a

    return c
    
def contrast_dark_light(color:Color):
    c:Color = Color(make_color(color))
    n = 1.5 if get_luminance(c) < 0.5 else 0.5
    for i in range(3):
        c[i] = int(min(max(c[i] * n, 0), 255))
    return c

def contrast_color(color:Color) -> Color:
    c:Color = Color(make_color(color))
    for i in range(3):
        c[i] = ((c[i] - int(255/3)) %256)

    return c

def max_bright(color:Color) -> Color:
    c:Color = Color(make_color(color))
    max_c = c[0]
    for i in range(1,3):
        if c[i] > max_c: max_c = c[i]
    for i in range(3):
        c[i] = int(c[i] * (255/max_c))
    return c

def step_to(from_color:Color, to_color:Color, luminance_only=False):
    fc:Color = Color(make_color(from_color))
    tc:Color = Color(make_color(to_color))
    if luminance_only:
        lum_f = get_luminance(fc)
        lum_t = get_luminance(tc)
        avg = (lum_f+lum_t)/2
        for i in range(4):
            fc[i] = int(fc[i]*avg)
    else:
        for i in range(4):
            fc[i] += int((tc[i]-fc[i])/2)
    return fc

#-- OTHER FUNCTIONS --#
def validate_color(color:tuple) -> Color | False:
    try:
        c = Color(color)
    except (ValueError, TypeError) as error:
        print(f"Error occured in Color.validate_color: {error}. color: {color}")
        return False
    return c

def make_color(color:Color|tuple) -> Color:
    c = color
    if callable(color):
        if isinstance(c(), Color): return c
        elif isinstance(c(), (tuple,list)) and len(c()) in [3,4]:
            return lambda: Color(*c())
        else:
            raise TypeError(f"Invalid color: {c}")
    else:
        if isinstance(c, Color): return c
        elif isinstance(c, (tuple,list)) and len(c) in [3,4]:
            return Color(*c)
        else:
            raise TypeError(f"Invalid color: {c}")

def get_luminance(color:Color) -> float:
    r,g,b,a = make_color(color).normalize()
    luminance = 0.2126*r + 0.7152*g + 0.0722*b
    return luminance


