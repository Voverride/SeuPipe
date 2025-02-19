from dataManager.users import *

def get_all_passcode()->list:
    """
    Retrieve all users' passcodes along with their permissions and statuses.

    This function collects the data of all users, including their passcode, permission status, 
    and account status. It organizes the data into a list of dictionaries, with options for 
    permissions and account statuses.

    Returns:
        list: A list of dictionaries containing user information such as ID, passcode, 
              permission settings, and account status.
    """
    users = get_users()
    data = [
        {
            'id': passcode.doc_id,
            'Passcode': passcode['usrname'],
            'Permission': {
                'options': [
                    {'label': 'Editable', 'value': 'Editable'},
                    {'label': 'ReadOnly', 'value': 'ReadOnly'},
                ],
                'allowClear': False,
                'bordered': False,
                'value': 'Editable' if passcode['usrtype'] > 0 else 'ReadOnly'
            },
            'Status': {
                'options': [
                    {'label': 'Active', 'value': 'Active'},
                    {'label': 'Disabled', 'value': 'Disabled'},
                ],
                'allowClear': False,
                'bordered': False,
                'value': 'Disabled' if passcode['disabled'] else 'Active'
            },
        } for passcode in users
    ]
    return data
def get_old_passcode(id:int)->str:
    """
    Retrieve the passcode (username) of a user by their ID.

    Args:
        id (int): The unique ID of the user whose passcode is being retrieved.

    Returns:
        str: The username (passcode) of the user.
    """
    return get_user_with_id(id)['usrname']

def update_passcode(id:int, new_passcode:str)->None:
    """
    Update the passcode (username) for a user by their ID.

    Args:
        id (int): The unique ID of the user whose passcode needs to be updated.
        new_passcode (str): The new passcode (username) to be set for the user.

    Returns:
        None
    """
    update_user_with_ids([id], new_name=new_passcode)

def update_permission_status(id:int, field:str, new_value:str)->None:
    """
    Update the permission or status field for a user based on the provided field and new value.

    Args:
        id (int): The unique ID of the user whose permission/status is being updated.
        field (str): The field to update. Can be 'Permission' or 'Status'.
        new_value (str): The new value to set for the field ('Editable', 'ReadOnly', 'Active', or 'Disabled').

    Returns:
        None
    """
    if field != 'Status':
        type = 0
        if new_value != 'ReadOnly':
            type = 1
        update_user_with_ids([id], new_type=type)
    else:
        disabled = True
        if new_value != 'Disabled':
            disabled = False
        update_user_with_ids([id], disabled=disabled)

def create_passcode(passcode:str, permission:str, status:str)->dict:
    """
    Create a new user with a specified passcode, permission, and status.

    This function checks if the passcode already exists, and if not, it adds the user with 
    the given details. It returns the user information as a dictionary, including the options 
    for permission and status fields.

    Args:
        passcode (str): The passcode (username) for the new user.
        permission (str): The permission level for the new user ('Editable' or 'ReadOnly').
        status (str): The status of the user ('Active' or 'Disabled').

    Returns:
        dict: A dictionary containing the newly created user's details, including options for 
              permission and status.
        None: If the passcode already exists, returns None.
    """
    users = search_user(usrname=passcode)
    if users:
        return None
    
    type = 0
    if permission == 'Editable':
        type = 1
    disabled = True
    if status == 'Active':
        disabled = False
    
    id = add_user(passcode, type, disabled)
    
    user_item = {
        'id': id,
        'Passcode': passcode,
        'Permission': {
            'options': [
                {'label': 'Editable', 'value': 'Editable'},
                {'label': 'ReadOnly', 'value': 'ReadOnly'},
            ],
            'allowClear': False,
            'bordered': False,
            'value': permission
        },
        'Status': {
            'options': [
                {'label': 'Active', 'value': 'Active'},
                {'label': 'Disabled', 'value': 'Disabled'},
            ],
            'allowClear': False,
            'bordered': False,
            'value': status
        },
    }
    return user_item

def remove_passcode(ids)->None:
    """
    Remove users from the system based on the provided user IDs.

    Args:
        ids (list): A list of user IDs to be removed.

    Returns:
        None
    """
    remove_user_with_ids(ids)