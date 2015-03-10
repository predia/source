class cstruct:
	def get(self,field):
		return self.__dict__[field]
	
	def __call__(self):
		print vars(self)
	
	def isfield(self,field):
		found = 0
		field_list = dir(self)
		for i in dir(self):
			if field == i:
				found = 1
		return found
		
	def fields(self):
		return vars(self).keys()
		
	def isSetTrue(self,field):
		found = 0
		field_list = dir(self)
		for i in dir(self):
			if field == i and self.get(i) == True:
				found = 1
		return found
		
	def get(self,field):
		return self.__dict__[field]	
	def __call__(self):
		print vars(self)
	
	def isfield(self,field):
		found = 0
		field_list = dir(self)
		for i in dir(self):
			if field == i:
				found = 1
		return found
		
	def get(self,field):
		return self.__dict__[field]
