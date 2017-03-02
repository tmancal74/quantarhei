# -*- coding: utf-8 -*-

match_number = r'(\d+(?:\.\d+)?)'
match_word = r'([^"]*)'

class FeatureFileGenerator:
    
    def __init__(self, feature):
        
        self._begin = "#begin_feature"
        self._end = "#end_feature"
        
        self._out = ""
        self._save = False
    
        self._feature = feature
        self._out += self._process_feature()
        self._has_feature = True
 
        self._examples = None
        self._has_examples = False
        
        
        self._givens = []
        self._has_givens = False
        self._whens = []
        self._has_whens = False
        self._thens = []
        self._has_thens = False
        
    def _process_feature(self):
        save_directly = True
        out = ""
        ex = ""
        for line in self._feature.splitlines():
            strip = line.strip()
            if strip == "Examples:":
                save_directly = False
            if save_directly:
                out += line+"\n"
            else:
                ex += line+"\n"
        self.examples = ex
        self._has_examples = True
        return out
            

    def _set_examples(self, ex):
        if self._has_thens:
            self._out += ex
        self._examples = ex
        self._has_examples = True        
        
    def _save_string(self, string):
        
        for line in string.splitlines():            
            strip = line.strip()
            if strip == self._begin:
                self._save = True
            elif strip == self._end:
                self._save = False
            else:
                if self._save:
                    self._out += line+"\n" 
                    
    def add_Given(self, func):
        if not self._has_feature:
            raise Exception()
        string = func.__doc__
        self._save_string(string)
        if self._save == True:
            raise Exception("Feature file content did not finish before end of the string")
        self._givens.append(func.__doc__)
        self._has_givens = True
        
    def add_When(self, func):
        if not self._has_givens:
            raise Exception()
        string = func.__doc__
        self._save_string(string)
        if self._save == True:
            raise Exception("Feature file content did not finish before end of the string")
        self._whens.append(func.__doc__)
        self._has_whens = True
        
    def add_Then(self, func):
        if not self._has_whens:
            raise Exception()
        string = func.__doc__
        self._save_string(string)
        if self._save == True:
            raise Exception("Feature file content did not finish before end of the string")
        self._thens.append(func.__doc__)
        self._has_thens = True
        
        
    def generate_feature_file(self, filename):
        self._set_examples(self.examples)
        #print(self._out)
        with open(filename,"w") as file:
            file.write(self._out)

