# -*- coding: utf-8 -*-

import tempfile
import os
import io
import tarfile
import shutil 

import numpy 
#import matplotlib.pyplot as plt
import xml.dom.minidom
#import scipy.interpolate
#from ..utils.types import Integer

class MatrixData:

    version = "0.1"
    date = "2015-04-09"


    def __init__(self,dims=(0),name="",data=""):
        
        self.name = name
        
        if data == "":
            self.data = numpy.zeros(dims)       
        else:
            self.data = data
        

    def setName(self,name):
        self.name = name

    # returns the name of the data set        
    def getName(self):
        return self.name
        
    def getDim(self,n):
        return self.data.shape[n]
        
    def getShape(self):
        return self.data.shape
        
    def getRank(self):
        return self.data.ndim
        
    def setData(self,data):
        self.data = data
        
        
        
        
    def saveBinaryData(self,filename):
        numpy.savez_compressed(filename,data=self.data)
        
    def loadBinaryData(self,filename):
        self.data = numpy.load(filename)["data"]
        
                
    def exportDataToText(self,filename):
        numpy.savetxt(filename,self.data)            
        
    def importDataFromText(self,filename):
        self.data = numpy.loadtxt(filename)
        
        
    def saveNoseFormat(self,name):

        filename = ("description.xml", "content.npz")

        ndir,nfn = os.path.split(name)

        with tempfile.TemporaryDirectory() as tdir:
        
            # description
            find = os.path.join(tdir,filename[0])
            fintd = open(find,"w")
            self.__writeDescription(fintd)
            fintd.close()

            find = os.path.join(tdir,filename[1]) 
            self.saveBinaryData(find)
            
            nfn = os.path.join(tdir,nfn)
            tf = tarfile.open(nfn,"w|gz")
            for fname in filename:
                lfname = os.path.join(tdir,fname)
                tf.add(lfname,arcname=fname)
            tf.close()    
        
            if not ndir == '':
                dest = os.path.join(os.getcwd(),ndir)
                if not os.path.exists(dest):
                    os.makedirs(dest)
                if not os.path.isdir(dest):
                    raise Exception
            else:
                dest = os.getcwd()
                
            shutil.move(tf.name,dest)    

    
    def loadNoseFormat(self,name):
        
        fd = open(name,"rb")

        fbytes = fd.read()
        file_like_object = io.BytesIO(fbytes)
        
        tar = tarfile.open(fileobj=file_like_object) 
        
        ii = 1
        for member in tar.getmembers():
            f = tar.extractfile(member)
            if ii == 1:
                cont = f.read().decode("utf-8")
                self.__readDescription(cont)
                ii += 1
            else:
                self.data = numpy.load(f)["data"]
                self.__compare_with_description()                
                
        fd.close()
    
    def __readDescription(self,desc):
        print(desc)
        shp = self.data.shape
        self.assert_shape(shp)
        tp = self.noseType
        self.assert_nose_type(tp)
        
        
    def __compare_with_description(self):
        pass
    
    
    def __writeDescription(self,fd,type_specific_data=None):
        
        dom = xml.dom.minidom.getDOMImplementation()
        tree = dom.createDocument(None,"description",None)
        root = tree.documentElement
        
        element = tree.createElement("nose-matrix-data")
        
        element.setAttribute("format-version",self.version)
        element.setAttribute("format-date",self.date)
        
        nameElement = tree.createElement("name")
        text_element = tree.createTextNode(self.name)
        nameElement.appendChild(text_element)
        element.appendChild(nameElement)
        
        typeElement = tree.createElement("nose-type")
        text_element = tree.createTextNode(self.noseType)
        typeElement.appendChild(text_element)
        element.appendChild(typeElement)

        shapeElement = tree.createElement("data-shape")
        text_element = tree.createTextNode(str(self.data.shape))
        shapeElement.appendChild(text_element)
        element.appendChild(shapeElement)
        
        fileElement = tree.createElement("data-file")
        text_element = tree.createTextNode("content.npz")
        fileElement.appendChild(text_element)
        element.appendChild(fileElement)        

        specElement = tree.createElement("type-specific-data")
        if type_specific_data == None:
            text_element = tree.createTextNode("None")
            specElement.appendChild(text_element)
        else:
            specElement.appendChild(type_specific_data)
            
        element.appendChild(specElement)         
        
        root.appendChild(element)        
        tree.writexml(fd, encoding="UTF-8",addindent="  ",newl="\n")
        
        
    def save(self,name):
        self.saveNoseFormat(name)
        
    def load(self,name):
        self.loadNoseFormat(name)

    def assert_shape(self,shp):
        print("Shape: Base class")
        
    def assert_nose_type(self,ntype):
        print("NType: Base class")
        
    def self_test(self):
        """

        Tests the functionality of the class
        
        """
 
        dd = MatrixData("data1")

        # generate test data


        dd.setData(numpy.ones((101,101)))
       
        print("Name: ",dd.getName())
        print("Rank: ",dd.getRank())
        print("First dim: ", dd.getDim(1))
        

        # test Txt export
        
        dd.exportDataToText("test.dat")
        
        cc = MatrixData("data2")
        cc.importDataFromText("test.dat")
        
        # compare original and saved data
        
        print("New dim: ",cc.getDim(1))


        # test saving with Nose Binary Data format

        fname = "data.nbd"
        try:
            cc.save(fname)
        except:
            os.remove(fname)
            cc.save(fname)
    
        aa = MatrixData("data3")
        aa.load(fname)

        print(aa.data)       





