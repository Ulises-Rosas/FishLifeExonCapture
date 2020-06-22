
def checkstep(step, step_choices):
    
    if step is None:
        sys.stdout.write("\nPlease, introduce an step\n")
        exit()
        
    if not step in step_choices:
        sys.stdout.write("\nPlease, choice between these steps:\n\n")
        sys.stdout.write("\n".join( [ "- " + i for i in step_choices] ) )
        exit()
