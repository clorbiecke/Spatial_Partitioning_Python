import pygame
from pygame import Rect
from pygame import Color, Surface, mouse, key
from pygame.math import Vector2
from pygame.event import Event
from physics_objects import PhysicsObject
import colors
import fonts

import itertools
import math


class Camera(Rect):
    def __init__(self, left:float, top:float, width:float, height:float, surface:Surface):
        super().__init__(left, top, width, height)
        self.surface:Surface = surface

    @property
    def pos(self) -> Vector2:
        return Vector2((self.left + self.width/2.0), (self.top + self.height/2.0))
    @pos.setter
    def pos(self, pos:Vector2):
        self.left = pos.x - self.width/2.0
        self.top = pos.y - self.height/2.0

    def center(self, center:Vector2):
        print(center.__class__.__name__)
        self.pos = Vector2(center)
        self.update(self.clamp(self.surface.get_rect()))
    
    def view(self) -> Rect:
        return Rect(self.left, self.top, self.width, self.height)
    
    def draw(self, window:Surface):
        x = (window.get_width()-self.width)*0.5
        y = (window.get_height()-self.height)*0.5
        window.blit(self.surface, (x,y), self)





class Bar():
    bar_anchors = {
        'left': (0, 0),
        'right': (1, 0),
        'top': (0, 0),
        'bot': (0, 1),
        'center': (0.5, 0.5)
    }
    fill_directions = {
        'vertical': (0, 1),
        'horizontal': (1, 0),
        'scale': (1, 1),
    }

    def __init__(self, pos: tuple, size: tuple, max_value: float=1.0, start_value: float=0, anchor_pos:str='left', fill_direction: str='horizontal', 
                 split: bool=False, border_color: tuple=(255,255,255), border_width: int=1, fill_color: tuple=(255,255,255)):
        self.pos = Vector2(pos)
        self.size = Vector2(size)
        self.max_value = max_value
        self.curr_value = start_value
        self.anchor_pos = Vector2(self.bar_anchors[anchor_pos])
        # ensure fill direction and anchor are compatible
        self._fill_rect = lambda: self.fill_rect_default()
        self._fill_direction = None
        self._split = None
        self.border_color = border_color
        self.border_width = border_width
        self.fill_color = fill_color        
        
        self.set_fill_direction(fill_direction)
        self.set_split(split)


    @property
    def fill_direction(self):
        return self._fill_direction
    def set_fill_direction(self, direction: str):
        self._fill_direction = Vector2(self.fill_directions[direction])

    @property
    def split(self):
        return self._split
    def set_split(self, split:bool=False):
        self._split = split
        if self.split == True:
            self._fill_rect = lambda: self.fill_rect_split()
        else:
            self._fill_rect = lambda: self.fill_rect_default()

    @property
    def fill_rect(self):
        return self._fill_rect()
    @property
    def fullness(self):
        n = self.curr_value/self.max_value
        if abs(n) > 1:
            n = n/abs(n)
        return n

    def set_curr_value(self, value: float):
        if abs(value) > self.max_value:
            self.curr_value = (value/abs(value))*self.max_value
        else:
            self.curr_value = value

    def draw(self, window):
        edge=Vector2(2,2)
        # draw border
        frame = pygame.Surface(self.size+2*edge, pygame.SRCALPHA)
        pygame.draw.rect(frame, self.border_color, frame.get_rect(), self.border_width)

        # frame.fill(self.border_color)
        # pygame.draw.rect(frame, self.border_color, (self.pos-edge, self.size+2*edge), self.border_width)
        # draw fill
        fill = pygame.Surface(self.fill_rect[1], pygame.SRCALPHA)
        pygame.draw.rect(fill, self.fill_color, fill.get_rect())

        # fill.fill(self.fill_color)
        # pygame.draw.rect(frame, self.fill_color, self.fill_rect)
        window.blit(frame, self.pos-edge)
        window.blit(fill, self.fill_rect[0])
    
    def rel_to_abs(self, vector):
        return vector.elementwise() * self.size.elementwise()
    
    def fill_rect_default(self):
        pos_c = (self.anchor_pos.elementwise() * self.fill_direction.elementwise()) + (-self.anchor_pos.elementwise() * self.fill_direction.elementwise()) * self.fullness# mult by fullness to get new fill pos
        size_c = -self.fill_direction * (1-self.fullness)# (-0.5, 0) * (1-0) = (-0.5, 0)
        pos = self.pos + self.rel_to_abs(pos_c)
        size = self.size + self.rel_to_abs(size_c) # (200, 40) + ((-0.5, 0)*(200, 40)) = (100, 0)
        return (pos, size)
    
    def fill_rect_split(self):
        fullness = self.fullness * (0.5 if self.fill_direction != Vector2(self.fill_directions['vertical']) else -0.5)
        anchor_pos = Vector2(0.5,0.5).elementwise() * self.fill_direction.elementwise()
        size_c = (-self.fill_direction * (1-abs(fullness)))
        pos_c = Vector2(0,0)

        if fullness < 0:
            pos_c = anchor_pos*fullness*2
        size = self.size + self.rel_to_abs(size_c)
        pos = self.pos + self.rel_to_abs(anchor_pos + pos_c)

        return (pos, size)
        pass


class Button():
    def __init__(self, txt:str, pos:Vector2, width:int, height:int, txt_color:Color|tuple, fill_color:Color|tuple, out_color:Color|tuple, on_click:lambda:(), on_hover:lambda:(), auto_txt_resize=True):
        self.txt = txt
        self.pos = pos # anchored to top left of rect
        
        self.width = width
        self.height = height
        self.fill_color = fill_color
        self.txt_color = txt_color
        self.out_color = out_color
        self.on_click = on_click
        self.on_hover = on_hover
        self.auto_txt_resize = auto_txt_resize

        self.txt_size = self.height/2
        self.txt_size_ratio = 0.6
        self.font = fonts.txt_size(self.txt_size)
        self.out_width = 3
        self.hovering = False
        self.pressed = False

    #- METHODS -#
    def update(self, dt):
        # print(f"Called update() on {self.txt} btn")
        # print(f"\tpos: {self.pos}, width: {self.width}, height: {self.height}, mouse_in_rect: {self.is_within_rect(mouse.get_pos())}, mpos: {mouse.get_pos()}")
        # TODO: add event poll to call on_click() and on_hover()
        if self.is_within_rect(mouse.get_pos()):
            self.hovering = True
            self.on_hover
            mdown = mouse.get_pressed()
            if mdown[0]:
                if not self.pressed:
                 self.pressed = True
                 self.on_click
        else:
            self.hovering = False
            self.pressed = False
        # pass

    def draw(self, surface:Surface):
        # create button surface
        fill_color = self.fill_color if not self.hovering else colors.contrast_dark_light(self.fill_color)
        out_color = self.out_color if not self.hovering else colors.contrast_dark_light(self.out_color)
        txt_color = self.txt_color if not self.hovering else colors.contrast_dark_light(self.txt_color)
        if self.pressed:
            fill_color = colors.dark(fill_color)
            out_color = colors.step_to(out_color, fill_color, True)
            txt_color = colors.step_to(txt_color, fill_color, True)

        button = Surface((self.width, self.height), pygame.SRCALPHA)
        button.fill(fill_color)

        # add button outline
        pygame.draw.rect(button, out_color, (0,0, button.get_width(), button.get_height()), self.out_width)

        # add button text
        txt = self.font.render(self.txt, True, txt_color)
        if self.auto_txt_resize: 
            c=self.txt_size
            n = self.resize_txt(txt.get_width(), txt.get_height())
            if c != n: txt = self.font.render(self.txt, True, txt_color)
        button.blit(txt, ((button.get_width()-txt.get_width())/2, (button.get_height()-txt.get_height())/2))

        # blit button to surface
        surface.blit(button, self.pos)
       
    def is_within_rect(self, pos:Vector2) -> bool:
        p = Vector2(pos)
        return (p.x > self.pos.x and p.x < self.pos.x+self.width and p.y > self.pos.y and p.y < self.pos.y+self.height)
    
    def resize_txt(self, width:int, height:int) -> int:
        w = (width/self.width)
        h = (height/self.height)
        r = self.txt_size_ratio
        if (max(w,h) > r) or (max(w,h) < r/2):
            factor = r / max(w,h)
            self.txt_size *= factor
            self.font = fonts.txt_size(self.txt_size)
        return self.txt_size 

    #- PROPERTIES -#
    @property
    def txt(self) -> str:
        return self._txt() if callable(self._txt) else self._txt
    @txt.setter
    def txt(self, txt:str):
        self._txt = txt

    @property
    def pos(self) -> Vector2:
        return Vector2(self._pos() if callable(self._pos) else self._pos)
    @pos.setter
    def pos(self, pos:Vector2):
        self._pos = pos

    @property
    def width(self) -> int:
        return int(self._width() if callable(self._width) else self._width)
    @width.setter
    def width(self, width:int):
        self._width = width
    
    @property
    def height(self) -> int:
        return int(self._height() if callable(self._height) else self._height)
    @height.setter
    def height(self, height:int):
        self._height = height

    @property
    def fill_color(self) -> Color:
        return self._fill_color() if callable(self._fill_color) else self._fill_color
    @fill_color.setter
    def fill_color(self, fill_color:Color|tuple):
        self._fill_color = fill_color

    @property
    def txt_color(self) -> Color:
        return self._txt_color() if callable(self._txt_color) else self._txt_color
    @txt_color.setter
    def txt_color(self, txt_color:Color|tuple):
        self._txt_color = txt_color

    @property
    def out_color(self) -> Color:
        return self._out_color() if callable(self._out_color) else self._out_color
    @out_color.setter
    def out_color(self, out_color:Color|tuple):
        self._out_color = out_color

    @property
    def on_click(self):
        self._on_click() if callable(self._on_click) else self._on_click
    @on_click.setter
    def on_click(self, on_click:lambda:()):
        self._on_click = on_click

    @property
    def on_hover(self):
        self._on_hover() if callable(self._on_hover) else self._on_hover
    @on_hover.setter
    def on_hover(self, on_hover:lambda:()):
        self._on_hover = on_hover

# TODO: add wrap text bool
class TextBox():
    horizontal_alignment = {
        "left" : 0,
        "center" : 0.5,
        "right" : 1
    }
    vertical_alignment = {
        "top" : 0,
        "center" : 0.5,
        "bottom" : 1
    }

    def __init__(self, txt:str, pos:Vector2, width:int, height:int, txt_color:Color|tuple, fill_color:Color|tuple, out_color:Color|tuple, 
                 auto_txt_resize=True, h_align="left", v_align="center"):
        self.txt = txt
        self.pos = pos # anchored to top left of rect
        
        self.width = width
        self.height = height
        self.fill_color = fill_color
        self.txt_color = txt_color
        self.out_color = out_color
        self.auto_txt_resize = auto_txt_resize
        self.h_align = h_align
        self.v_align = v_align

        self.txt_size = self.height/2
        self.txt_size_ratio = 0.6
        self.font = fonts.txt_size(self.txt_size)
        self.out_width = 3

    #- METHODS -#
    def update(self, dt=0):
        pass

    def draw(self, surface:Surface):
        # create button surface
        fill_color = self.fill_color# if not self.hovering else colors.contrast_dark_light(self.fill_color)
        out_color = self.out_color# if not self.hovering else colors.contrast_dark_light(self.out_color)
        txt_color = self.txt_color# if not self.hovering else colors.contrast_dark_light(self.txt_color)
        # if self.pressed:
        #     fill_color = colors.dark(fill_color)
        #     out_color = colors.step_to(out_color, fill_color, True)
        #     txt_color = colors.step_to(txt_color, fill_color, True)
        box = Surface((self.width, self.height), pygame.SRCALPHA)
        box.fill(fill_color)

        # add text box outline
        pygame.draw.rect(box, out_color, (0,0, box.get_width(), box.get_height()), self.out_width)

        # add text box text
        txt = self.font.render(self.txt, True, txt_color)
        if self.auto_txt_resize: 
            c=self.txt_size
            n = self.resize_txt(txt.get_width(), txt.get_height())
            if c != n: txt = self.font.render(self.txt, True, txt_color)
        
        pos = Vector2(0)
        pos.x = self.h_align * (box.get_width()-txt.get_width())
        pos.y = self.v_align * (box.get_height()-txt.get_height())
        box.blit(txt, pos)

        # blit button to surface
        surface.blit(box, self.pos)
       
    def is_within_rect(self, pos:Vector2) -> bool:
        p = Vector2(pos)
        return (p.x > self.pos.x and p.x < self.pos.x+self.width and p.y > self.pos.y and p.y < self.pos.y+self.height)
    
    def resize_txt(self, width:int, height:int) -> int:
        w = (width/self.width)
        h = (height/self.height)
        r = self.txt_size_ratio
        if (w > r or h > r):
            factor = r / max(w,h)
            self.txt_size *= factor
            self.font = fonts.txt_size(self.txt_size)
        elif(w < r/2 or h < r/2):
            factor = r / max(w,h)
            self.txt_size *= factor
            self.font = fonts.txt_size(self.txt_size)
        return self.txt_size 

    
        
    #- PROPERTIES -#
    @property
    def txt(self) -> str:
        s = self._txt() if callable(self._txt) else self._txt
        if isinstance(s, str): return s
        else: return f"{s}"
         
    @txt.setter
    def txt(self, txt:str):
        self._txt = txt


    @property
    def pos(self) -> Vector2:
        return Vector2(self._pos() if callable(self._pos) else self._pos)
    @pos.setter
    def pos(self, pos:Vector2):
        self._pos = pos

    @property
    def width(self) -> int:
        return int(self._width() if callable(self._width) else self._width)
    @width.setter
    def width(self, width:int):
        self._width = width
    
    @property
    def height(self) -> int:
        return int(self._height() if callable(self._height) else self._height)
    @height.setter
    def height(self, height:int):
        self._height = height

    @property
    def fill_color(self) -> Color:
        return self._fill_color() if callable(self._fill_color) else self._fill_color
    @fill_color.setter
    def fill_color(self, fill_color:Color|tuple):
        self._fill_color = fill_color

    @property
    def txt_color(self) -> Color:
        return self._txt_color() if callable(self._txt_color) else self._txt_color
    @txt_color.setter
    def txt_color(self, txt_color:Color|tuple):
        self._txt_color = txt_color

    @property
    def out_color(self) -> Color:
        return self._out_color() if callable(self._out_color) else self._out_color
    @out_color.setter
    def out_color(self, out_color:Color|tuple):
        self._out_color = out_color

    @property
    def txt_size(self) -> int:
        return self._txt_size
    #TODO: make txt size setter to ensure it is an int
    @txt_size.setter
    def txt_size(self, size:int):
        self._txt_size = int(size)

    @property
    def h_align(self) -> float:
        return self._h_align
    @h_align.setter
    def h_align(self, h_align:str):
        if h_align in self.horizontal_alignment:
            self._h_align = self.horizontal_alignment.get(h_align)
        else:
            raise ValueError(f"Error: h_align must be a value in {self.horizontal_alignment.keys()}. Given value: {h_align}.")

    @property
    def v_align(self) -> str:
        return self._v_align
    @v_align.setter
    def v_align(self, v_align:str):
        if v_align in self.vertical_alignment:
            self._v_align = self.vertical_alignment.get(v_align)
        else:
            raise ValueError(f"Error: v_align must be a value in {self.vertical_alignment.keys()}. Given value: {v_align}.")






