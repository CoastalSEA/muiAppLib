Installing the XXXX App
XXXX is installed as an App and requires muitoolbox and dstoolbox to be installed. The download for each of these includes the code, documentation and example files. The files required are: 
> dstoolbox.mltbx
> muitoolbox.mltbx
> XXXX.mlappinstall

Installing the Toolboxes

You can check whether the required toolboxes are already installed using:
addons = matlab.addons.installedAddons
which displays the installation details in the command window (you may
need to scroll the table to see both toolboxes)

If not already installed, the two toolboxes can be installed using the Add-Ons>Manage Add-Ons option on the HOME tab of Matlab™. Alternatively, right-click the mouse on the ‘mltbx’ files and select Install. All the folder paths are initialised upon installation and the location of the code is also handled by Matlab™. The location of the code can be found using the options in the Manage Add-Ons UI.

Installing the App

The App is installed using the Install Apps button on the APPS tab in Matlab™. Alternatively, right-click the mouse on the ‘mlappinstall’ file and select Install. Again all the folder paths are initialised upon installation and the location of the code is handled by Matlab™.

Once installed, the App can be run from the APPS tab. This sets the App environment paths, after which the App can be run from the Command Window using:
>>  XXXX;

The App environment paths can be saved using the Set Path option on the Matlab™ HOME tab.
Documentation can be viewed from the Supplementary Software in the Matlab™ Documentation. The location of the code can be accessed by hovering over the App icon and then finding the link in the pop-up window.

The following line of code will install the XXXX App
>>a ppinfo = matlab.apputil.install('XXXX.mlappinstall')     
which displays the installation details in the command window

Once the App is installed, to check what Apps are installed, use:
>> appinfo = matlab.apputil.getInstalledAppInfo