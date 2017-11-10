# -*- coding: utf-8 -*-
"""
    IPython magic commands and their equivalents for Python console
    
    
    Magic commands
    --------------
    
    %template {list, fetch} [template_name]
        Lists available templates, or fetches a template by its name
        
        
        
    Console functions
    -----------------
    
    fetch_template(template_name)
        Fetches a template by its name
        
    list_templates()
        Lists available templates
        


"""
import pkg_resources


try:
    
    # if IPython is not present, 
    # we fail here and continue in the except section
    ip = get_ipython()
    _have_ip_ = True
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
            
            list_templates()
            
        elif command == "fetch":
            
            try:
                tname = words[1]
            except IndexError:
                print("Error: ``fetch`` command requires an argument\n")
                print("*** Printing __doc__ string ***\n")
                _print_template_help()
                return
            
            fetch_template(tname)

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

    
    
except NameError:
    _have_ip_ = False

    
    
def _print_template_help():
    """Prints doc string
    
    """
    print(template.__doc__)


def list_templates():
    print("""
List of templates
-----------------
    
ExcitonDynamics_Lindblad :
    Exciton dynamics calculated with relaxation specified by Lindblad form 
    
    """)

    
    
def fetch_template(tname):
    """Reads template file from the template directory
    
    
    """
    templates = {}
    templates["ExcitonDynamics_Lindblad"] = "excitondynamics_lindblad.py"
    
    resource_package = "quantarhei"  
    try:
        tfile = templates[tname]
    except KeyError:
        raise Exception("Template error: "+tname+" no such template")
        
    resource_path = '/'.join(('wizard', 'templates', tfile))  

    template = pkg_resources.resource_string(resource_package, resource_path)
     
    if _have_ip_:
        
        # if we have IPython, we put the template into the next cell
        ip = get_ipython()
        ip.set_next_input(template.decode("utf-8"))
        
    else:
        print(template.decode("utf-8"))

    
          
        


