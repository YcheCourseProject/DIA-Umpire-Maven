import os

if __name__ == '__main__':
    mount_dir_path = os.path.abspath(os.pardir)
    os.system('docker run  -v ' + mount_dir_path + ':/opt -it yche/yche-biocontainer /bin/bash -c "cd /opt/; zsh"')
