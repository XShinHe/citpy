#!/usr/bin/env python3
# coding=utf-8


class Atom:
	def __init__(self, **args):
		if 'name' in args:
			self.name = args['name']
		else:
			self.name = 'undef'
		if 'index' in args:
			self.index = args['index']
		else:
			self.index = 'undef'
		self.__autocomplete()
		
	def __autocomplete(self):
		self.mass = 1
		self.zval = 2
		self.property = 3
	
	def __repr__(self):
		return '<Atom>\nname  = %s\nindex = %s'%(self.name, self.index)
	def __str__(self):
		return 'name  = %s\nindex = %s'%(self.name, self.index)
		
	
class Mole:
	def __init__(self, **args):
		if 'name' in args:
			self.name = args['name']
		if 'formula' in args:
			self.formula = args['formula']
		if 'alist' in args:
			self.alist = args['alist']
			
	def __autotopo(self):
		


if __name__ == '__main__':
	a = Atom(name='He')
	print(a)
		
