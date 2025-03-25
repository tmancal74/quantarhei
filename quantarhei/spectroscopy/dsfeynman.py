# -*- coding: utf-8 -*-

class DSFeynmanDiagram():
    """ Double-sided Feynman diagrams


    """
    
    def __init__(self, ptype="non-defined"):

        self.states = {}
        self.pointer = 0
        self._pic_rep = "\n"
        self.finished = False

        self.type = ptype
        
        self.states[self.pointer] = ["g", "g"]
        self.add_states_line()


    def add_states_line(self):
        self._pic_rep = "      | "+self.states[self.pointer][0]+ \
                        "      "+self.states[self.pointer][1]+" |\n"+self._pic_rep


    def add_arrow(self, side, dir, to="a"):
        """Adds an arrow in a given direction and transits to a new state

        """
        self.pointer += 1
        if side == "left":
            self.states[self.pointer] = [to, self.states[self.pointer-1][1]]
            if dir == "--->":
                self._pic_rep = "  --->|----------|\n"+self._pic_rep
            elif dir == "<---":
                self._pic_rep = "  <---|----------|\n"+self._pic_rep
                
        elif side == "right":
            self.states[self.pointer] = [self.states[self.pointer-1][0],to]
            if dir == "<---":
                self._pic_rep = "      |----------|<---\n"+self._pic_rep
            elif dir == "--->":
                self._pic_rep = "      |----------|--->\n"+self._pic_rep
                
        self.add_states_line()


    def finish(self, end="g"):
        """Finishes the diagram with a given state

        """
        self.add_arrow("left", "<---", end)
        self._pic_rep = "\n"+self._pic_rep
        self.count = self.pointer - 1
        self.finished = True


    def _check_finished(self):
        """Checkes if the diagram is finished

        """
        if not self.finished:
            raise Exception("An unfinisged diagram: finish() method has to be called first")


    def evolution_operators(self):

        self._check_finished()
            
        kk = 0
        Uops = ""
 
        for key in self.states:
            if kk > 0 and kk < self.count+1:
                sts = self.states[key]
                rightstate = sts[1]
                if rightstate == "g":
                    Uops += "Ugd(t"+str(kk)+")"
                else:
                    Uops += "Ued("+rightstate+",t"+str(kk)+")"
                if kk < self.count:
                    Uops += "*"
            kk += 1

        kk = 0
        rUop = ""
        
        for key in self.states:
             if kk > 0 and kk < self.count+1:
                sts = self.states[key]
                leftstate = sts[0]
                if leftstate == "g":
                    rUop = "Ug(t"+str(kk)+")"+rUop
                else:
                    rUop = "Ue("+leftstate+",t"+str(kk)+")"+rUop
                if kk < self.count:
                    rUop = "*"+rUop
             kk += 1

        Uops += "*"+rUop
        
        return Uops
            
    
    def __str__(self):
        self._check_finished()
        return self._pic_rep

    
    def report(self):
        self._check_finished()
        print("\nDiagram of",self.type,"type\n")
        print(self)
        print("Light interaction count:", self.count)   


class R1g_Diagram(DSFeynmanDiagram):
    """R1g diagram

    Diagram of R1g type
    
    
          | g      g |
      <---|----------|
          | a      g |
          |----------|--->
          | a      b |
          |----------|<---
          | a      g |
      --->|----------|
          | g      g |

    """
    
    def __init__(self, states=["a","b"]):
        super().__init__(ptype="R1g")
        self.add_arrow("left", "--->", to=states[0])
        self.add_arrow("right", "<---", to=states[1])
        self.add_arrow("right", "--->", "g" )
        self.finish()



class R2g_Diagram(DSFeynmanDiagram):
    """R2g diagram

    Diagram of R2g type
    
    
          | g      g |
      <---|----------|
          | b      g |
          |----------|--->
          | b      a |
      --->|----------|
          | g      a |
          |----------|<---
          | g      g |

    """

    def __init__(self, states=["a","b"]):
        super().__init__(ptype="R2g")
        self.add_arrow("right", "<---", to=states[0])
        self.add_arrow("left", "--->", to=states[1])
        self.add_arrow("right", "--->", "g" )
        self.finish()


class R3g_Diagram(DSFeynmanDiagram):
    """R3g diagram

    Diagram of R3g type
    
    
          | g      g |
      <---|----------|
          | b      g |
      --->|----------|
          | g      g |
          |----------|--->
          | g      a |
          |----------|<---
          | g      g |

    """

    def __init__(self, states=["a","b"]):
        super().__init__(ptype="R3g")
        self.add_arrow("right", "<---", to=states[0])
        self.add_arrow("right", "--->", to="g")
        self.add_arrow("left", "--->", to=states[1])
        self.finish()
    

class R4g_Diagram(DSFeynmanDiagram):
    """R4g diagram

    Diagram of R4g type
    
    
          | g      g |
      <---|----------|
          | b      g |
      --->|----------|
          | g      g |
      <---|----------|
          | a      g |
      --->|----------|
          | g      g |

    """

    def __init__(self, states=["a","b"]):
        super().__init__(ptype="R4g")
        self.add_arrow("left", "--->", to=states[0])
        self.add_arrow("left", "<---", to="g")
        self.add_arrow("left", "--->", to=states[1])
        self.finish()


class R1f_Diagram(DSFeynmanDiagram):
    """R1f diagram

    Diagram of R1f type
    
    
          | b      b |
      <---|----------|
          | f      b |
      --->|----------|
          | a      b |
          |----------|<---
          | a      g |
      --->|----------|
          | g      g |

    """

    def __init__(self, states=["a","b","f"]):
        super().__init__(ptype="R1f")
        self.add_arrow("left", "--->", to=states[0])
        self.add_arrow("right","<---", to=states[1])
        self.add_arrow("left", "--->", to=states[2])
        self.finish(end="b")
        

class R2f_Diagram(DSFeynmanDiagram):
    """R2f diagram

    Diagram of R2f type
    
          | a      a |
      <---|----------|
          | f      a |
      --->|----------|
          | b      a |
      --->|----------|
          | g      a |
          |----------|<---
          | g      g |
          
    """

    def __init__(self, states=["a","b","f"]):
        super().__init__(ptype="R2f")
        self.add_arrow("right", "<---", to=states[0])
        self.add_arrow("left","--->", to=states[1])
        self.add_arrow("left", "--->", to=states[2])
        self.finish(end="a")
