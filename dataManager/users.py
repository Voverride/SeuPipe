from tinydb import TinyDB
from tinydb import Query
from .alignment_d import alidata
import getpass
import os

directory = 'database/'
if not os.path.exists(directory):
    os.makedirs(directory)
db = TinyDB('database/db.json')
table = db.table('users')

admin = getpass.getuser()

q = Query()

def get_users()->list:
    """
    Retrieve a list of users from the database who are not the current admin.

    Returns:
        list: A list of users that do not have the same username as the admin.
    """
    return table.search(q.usrname != admin)

def add_user(usrname:str, usrtype:int, disabled:bool=False, usrhost:str=False)->int:
    """
    Add a new user to the 'users' table in the database.

    Args:
        usrname (str): The username(passcode) of the new user.
        usrtype (int): The user permission (e.g., 0->ReadOnly, 1->Editable).
        disabled (bool, optional): Indicates whether the user is disabled. Defaults to False.
        usrhost (str, optional): The host associated with the user. Defaults to False.

    Returns:
        int: The ID of the newly inserted user.
    """
    id = table.insert({'usrname': usrname, 'usrtype': usrtype, 'usrhost': usrhost, 'disabled': disabled})
    return id

def search_user(usrname:str=None, usrhost:str=None)->list:
    """
    Search for users based on the provided username and/or host.

    Args:
        usrname (str, optional): The username(passcode) of the user to search for. Defaults to None.
        usrhost (str, optional): The host of the user to search for. Defaults to None.

    Returns:
        list: A list of users that match the search criteria.
    """
    if usrname and usrhost:
        return table.search(q.usrname == usrname and q.usrhost == usrhost)
    elif usrname:
        return table.search(q.usrname == usrname)
    elif usrhost:
        return table.search(q.usrhost == usrhost)
    return []

def get_user_with_id(id:int)->dict:
    """
    Retrieve a user from the database by their unique ID.

    Args:
        id (int): The ID of the user to retrieve.

    Returns:
        dict: The user record with the specified ID.
    """
    return table.get(doc_id=id)

def remove_user_with_ids(ids:list)->None:
    """
    Remove users from the database using their IDs.

    Args:
        ids (list): A list of IDs corresponding to the users to be removed.
    """
    users = []
    for id in ids:
        users.append(get_user_with_id(id))
    table.remove(doc_ids=ids)
    for user in users:
        usrname = user['usrname']
        alidata.delete_user_data(usrname)

def update_user_with_ids(ids:list, new_name:str=None, new_type:str=None, new_host:str=None, disabled:bool=None)->None:
    """
    Update the user data for a list of user IDs with new values.

    Args:
        ids (list): The list of user IDs to update.
        new_name (str, optional): The new username(passcode) to set. Defaults to None.
        new_type (str, optional): The new user type to set. Defaults to None.
        new_host (str, optional): The new host to set. Defaults to None.
        disabled (bool, optional): The new disabled status. Defaults to None.
    """
    new_data = {}
    if new_name is not None:
        new_data['usrname'] = new_name
    if new_type is not None:
        new_data['usrtype'] = new_type
    if new_host is not None:
        new_data['usrhost'] = new_host
    if disabled is not None:
        new_data['disabled'] = disabled
    table.update(new_data, doc_ids=ids)