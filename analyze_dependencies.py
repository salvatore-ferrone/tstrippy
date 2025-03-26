import os
import re
import ast
from collections import defaultdict
import sys 
def find_imports(directory):
    imports = defaultdict(set)
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.py'):
                filepath = os.path.join(root, file)
                try:
                    with open(filepath, 'r') as f:
                        content = f.read()
                    
                    # Parse the file
                    tree = ast.parse(content)
                    
                    # Extract imports
                    for node in ast.walk(tree):
                        if isinstance(node, ast.Import):
                            for name in node.names:
                                imports[name.name].add(filepath)
                        elif isinstance(node, ast.ImportFrom):
                            if node.module is not None:
                                base_module = node.module.split('.')[0]
                                imports[base_module].add(filepath)
                except:
                    print(f"Could not parse {filepath}")
    
    return imports

if __name__ == "__main__":
    imports = find_imports("./tstrippy")
    
    # Filter out standard library and local imports
    standard_libs = set(sys.builtin_module_names)
    standard_libs.update(['os', 'sys', 're', 'math', 'datetime', 'time', 'random', 'json', 'csv'])
    
    print("External dependencies:")
    for module, files in sorted(imports.items()):
        # Skip standard library
        if module in standard_libs or module.startswith('tstrippy'):
            continue
        
        print(f"- {module}: used in {len(files)} files")