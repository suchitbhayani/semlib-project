import os 
import shutil
'''
Utility functions for project.
'''

def clean_files():
    '''
    Deletes all built files. Reverts to clean repo.
    '''
    data_dir = "data"
    if not os.path.exists(data_dir):
        print("`data/` directory does not exist.")
        return

    for root, dirs, files in os.walk(data_dir):
        for file in files:
            file_path = os.path.join(root, file)
            os.remove(file_path)

        for dir_name in dirs:
            dir_path = os.path.join(root, dir_name)
            shutil.rmtree(dir_path)

    print("All files and folders inside `data/` have been deleted.")