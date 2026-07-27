# GlennHT-GUI
A web-based graphical user interface to interact with plot3d files.

# Requirements
---
**Note about Ubuntu**

If you are using ubuntu make sure you first install build essentials. https://linuxize.com/post/how-to-install-gcc-on-ubuntu-20-04/ 
Some of the packages require gcc.

The front end requires node js. I reccomend following this guide https://docs.microsoft.com/en-us/windows/dev-environment/javascript/nodejs-on-wsl

---

- Tested on [Python 3.9.5](https://www.python.org/downloads/)
- Install packages with pip
    - Navigate to the backend folder: ```cd gui/backend/```
    - ```pip install -r requirements.txt```
- Tested on [Node v14.16.0](https://nodejs.org/en/download/)
- Install node modules
    - Navigate to the frontend folder: ```cd gui/frontend/```
    - ```npm install```

# Run backend
- Navigate to the backend folder: ```cd gui/backend/```
- ```python main.py```

# Run frontend
- Navigate to the frontend folder: ```cd gui/frontend/```
- ```npm i```
- ```npm start```

# View GUI
- Go to ```http://localhost:3000/``` in a web browser

# Interact with GUI
- Download mesh file: ```wget https://nasa-public-data.s3.amazonaws.com/plot3d_utilities/PahtCascade-ASCII.xyz```
- Upload file and view mesh via GUI.

# View API docs
- Go to ```http://localhost:5000/docs``` in a web browser
