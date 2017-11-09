# -*- coding: utf-8 -*-

import pkg_resources

from IPython.core.magic import register_line_magic

@register_line_magic
def template(line):
    """Quantarhei template magic command
    
    Commands
    ---------
    
    list
        Lists the available templates
        
    fetch template_name
        Opens a new cell and writes a template of given template_name
        into it
    
    
    """
    def _print_template_help():
        print(template.__doc__)
    
    def _print_template_list():
        print("""
List of templates
-----------------
        
ExcitonDynamics_Lindblad :
    Exciton dynamics calculated with relaxation specified by Lindblad form 
        
        """)
    
    def _fetch_template(tname):
        """Reads template file from the template directory
        
        """
        templates = {}
        templates["ExcitonDynamics_Lindblad"] = "excitondynamics_lindblad.py"
        
        resource_package = "quantarhei"  
        try:
            tfile = templates[tname]
        except KeyError:
            raise Exception("Template error: "+tname+" no such template")
            
        resource_path = '/'.join(('templates', tfile))  

        template = pkg_resources.resource_string(resource_package, resource_path)
            
        ip = get_ipython()
        ip.set_next_input(template.decode("utf-8"))

    
    # split the line of the magic command
    words = line.split()
    
    # get command
    try:
        command = words[0]
    except IndexError:
        print("Error: No command specified")
        print("Usage:")
        print("%template command [args]\n")
        print("*** Printing __doc__ string ***\n")
        _print_template_help()
        return
    
    # process the command
    if command == "list":
        
        _print_template_list()
        
    elif command == "fetch":
        
        try:
            tname = words[1]
        except IndexError:
            print("Error: ``fetch`` command requires an argument\n")
            print("*** Printing __doc__ string ***\n")
            _print_template_help()
            return
        
        _fetch_template(tname)
        
        

@register_line_magic
def example(line):
    """Quantarhei example magic command
    
    Commands
    ---------
    
    list
        Lists the available examples
        
    fetch example_name
        Opens a new cell and writes an example of given example_name
        into it
    
    
    """    
    words = line.split()
    for w in words:
        print(w)    


