import pam
from dataManager.users import *
from dash import set_props
from flask import request
from .notice import set_head_notice


def verify_modify_permission()->bool:
    """
        Verifies whether the current user has permission to modify.

        The function checks if the current user (based on their host) is valid, not disabled,
        and has an appropriate user type (greater than 0). If these conditions are met, 
        permission to modify is granted. Otherwise, a warning message is set and permission is denied.

        Returns:
            bool: True if the user has modify permission, False otherwise.
    """
    host = get_host()
    user = search_user(usrhost=host)
    if user and not user[0]['disabled'] and user[0]['usrtype']>0:
        return True
    else:
        set_head_notice('Permission Denied.', 'warning')
        return False

def get_host()->None:
    """
    Retrieve the remote address of the client making the request.

    Returns:
        str: The IP address of the client (remote host).
    """
    return request.remote_addr

def verify_user(usrname:str, password:str=None)->bool:
    """
    Verifies the user credentials and manages user access and session properties.

    This function checks whether the given username corresponds to an admin or a regular user,
    authenticates using the PAM (Pluggable Authentication Module) for the admin, and validates
    the status of regular users in the system. Upon successful verification, the session properties
    are updated, including the user ID, username, and the host from which the user is connecting.

    Args:
        usrname (str): The username/passcode of the user attempting to log in.
        password (str, optional): The password provided for authentication. Defaults to None.

    Returns:
        bool: True if the user is successfully authenticated, False otherwise.
    """
    host = get_host()
    access = False
    user = search_user(usrname=usrname)
    
    if usrname == admin:
        p = pam.pam()
        pwd = '' if password is None else password
        access = p.authenticate(usrname, pwd)
    else:
        if user and not user[0]['disabled']:
            access = True
    
    if access:
        userid = user[0].doc_id
        set_props('userid', dict(data=userid))
        set_props('Spatpy', dict(style=None))
        set_props('login-box', dict(visible=False))
        set_props('main-title-username', dict(children=usrname))
        update_user_with_ids([userid], new_host=host)
    
    return access

def verify_host()->None:
    """
    Verifies the host and updates the session properties for the user.

    This function checks whether the admin user exists in the database and creates one if not.
    It also checks the current host and determines whether the user should be allowed to log in
    based on their status and host. The session properties are updated accordingly.

    If the user is not found or disabled, the login box is shown, otherwise, the user's information
    is loaded into the session.

    Returns:
        None
    """
    if not search_user(usrname=admin):
        add_user(admin, 2)
    
    host = get_host()
    user = search_user(usrhost=host)
    
    if not user or user[0]['disabled']:
        set_props('login-box', dict(visible=True))

def restore_usrinfo()->None:
    """
    还原用户身份
    """
    host = get_host()
    user = search_user(usrhost=host)
    if user and not user[0]['disabled']:
        usrname = user[0]['usrname']
        userid = user[0].doc_id
        set_props('userid', dict(data=userid))
        set_props('main-title-username', dict(children=usrname))
        set_props('Spatpy', dict(style=None))

def logout(id:int)->None:
    """
    Logs out a user by removing the host information from their user record.

    Args:
        id (int): The ID of the user to log out.

    Returns:
        None
    """
    update_user_with_ids([id], new_host=False)