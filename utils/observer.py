# from typing import Dict, Callable, Any
# from collections.abc import MutableMapping

# class DictObserver:
#     """
#     字典变化观察器 - 单例模式
#     """
#     _instance = None
#     _initialized = False
    
#     def __new__(cls):
#         if cls._instance is None:
#             cls._instance = super(DictObserver, cls).__new__(cls)
#         return cls._instance
    
#     def __init__(self):
#         if not self._initialized:
#             self._watched_dicts = {}  # {id(original_dict): info}
#             self._wrapped_to_original = {}  # {id(wrapped_dict): id(original_dict)}
#             self._initialized = True
    
#     def observe(self, target_dict: Dict, callback: Callable[[str, Any, Any], None]):
#         """
#         监听字典变化
        
#         Args:
#             target_dict: 要监听的目标字典
#             callback: 回调函数，接收三个参数：操作类型('set', 'delete', 'update')、键、值
#         """
#         dict_id = id(target_dict)
        
#         # 如果字典还没有被包装过，创建一个代理对象
#         if dict_id not in self._watched_dicts:
#             wrapped_dict = ObservableDict(target_dict, self._notify_callbacks)
#             wrapped_id = id(wrapped_dict)
            
#             self._watched_dicts[dict_id] = {
#                 'original': target_dict,
#                 'wrapped': wrapped_dict,
#                 'callbacks': set()
#             }
#             self._wrapped_to_original[wrapped_id] = dict_id
        
#         # 添加回调函数
#         self._watched_dicts[dict_id]['callbacks'].add(callback)
        
#         # 返回可观察的代理字典
#         return self._watched_dicts[dict_id]['wrapped']
    
#     def disobserve(self, target_dict: Dict, callback: Callable = None):
#         """
#         解除监听
        
#         Args:
#             target_dict: 要解除监听的字典（可以是原始字典或包装后的字典）
#             callback: 要移除的特定回调函数，如果为None则移除所有回调
#         """
#         # 判断传入的是不是包装后的字典
#         dict_id = None
#         if isinstance(target_dict, ObservableDict):
#             wrapped_id = id(target_dict)
#             if wrapped_id in self._wrapped_to_original:
#                 dict_id = self._wrapped_to_original[wrapped_id]
#         else:
#             # 传入的是原始字典
#             dict_id = id(target_dict)
        
#         if dict_id is not None and dict_id in self._watched_dicts:
#             watched_info = self._watched_dicts[dict_id]
            
#             if callback is None:
#                 # 移除所有回调
#                 watched_info['callbacks'].clear()
#             else:
#                 # 移除特定回调
#                 watched_info['callbacks'].discard(callback)
            
#             # 如果没有回调了，清理整个监听
#             if len(watched_info['callbacks']) == 0:
#                 # 清理反向映射
#                 wrapped_id = id(watched_info['wrapped'])
#                 if wrapped_id in self._wrapped_to_original:
#                     del self._wrapped_to_original[wrapped_id]
#                 del self._watched_dicts[dict_id]
    
#     def _notify_callbacks(self, target_dict: Dict, operation: str, key: Any, value: Any):
#         """通知所有回调函数"""
#         dict_id = id(target_dict)
        
#         if dict_id in self._watched_dicts:
#             callbacks = self._watched_dicts[dict_id]['callbacks'].copy()
#             for callback in callbacks:
#                 try:
#                     callback(operation, key, value)
#                 except Exception as e:
#                     print(f"Callback error: {e}")


# class ObservableDict(MutableMapping):
#     """
#     可观察的字典包装器
#     """
    
#     def __init__(self, original_dict: Dict, notify_func: Callable):
#         self._data = original_dict
#         self._notify = notify_func
#         self._original_dict = original_dict
    
#     def __getitem__(self, key):
#         return self._data[key]
    
#     def __setitem__(self, key, value):
#         old_value = self._data.get(key, None) if key in self._data else None
#         is_new = key not in self._data
#         self._data[key] = value
        
#         if is_new:
#             self._notify(self._original_dict, 'set', key, value)
#         else:
#             self._notify(self._original_dict, 'update', key, (old_value, value))
    
#     def __delitem__(self, key):
#         value = self._data[key]
#         del self._data[key]
#         self._notify(self._original_dict, 'delete', key, value)
    
#     def __iter__(self):
#         return iter(self._data)
    
#     def __len__(self):
#         return len(self._data)
    
#     def __repr__(self):
#         return f"ObservableDict({self._data})"
    
#     def update(self, *args, **kwargs):
#         """
#         重写update方法以支持批量更新通知
#         """
#         updates = {}
        
#         if args:
#             other = args[0]
#             if hasattr(other, "items"):
#                 for key, value in other.items():
#                     updates[key] = value
#             else:
#                 for key, value in other:
#                     updates[key] = value
        
#         for key, value in kwargs.items():
#             updates[key] = value
        
#         for key, value in updates.items():
#             self.__setitem__(key, value)

# observer = DictObserver()

# if __name__ == "__main__":
#     """
#     测试代码
#     """
#     test_dict = {"name": "Alice", "age": {'a':25, 'b':30}}
#     def my_callback(operation, key, value):
#         print(f"字典发生变化: {operation}, 键: {key}, 值: {value}")
        
#     observed_dict = observer.observe(test_dict, my_callback)
        
#     print("初始字典:", observed_dict)
        
#     observed_dict["city"] = "Beijing"
        
#     observed_dict["age"] = 26
        
#     del observed_dict["name"]

#     observed_dict.update({"country": "China", "job": "Engineer"})
        
#     print("最终字典:", observed_dict)
        
#     observer.disobserve(observed_dict, my_callback)
        
#     observed_dict["test"] = "no callback"


from typing import Dict, Callable, Any, Union, Optional
from collections.abc import MutableMapping
import asyncio
import threading


class DictObserver:
    """
    字典变化观察器 - 单例模式（支持异步回调）
    """
    _instance = None
    _initialized = False
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(DictObserver, cls).__new__(cls)
        return cls._instance
    
    def __init__(self):
        if not self._initialized:
            self._watched_dicts = {}  # {id(original_dict): info}
            self._wrapped_to_original = {}  # {id(wrapped_dict): id(original_dict)}
            self._loop = None  # 存储主线程的事件循环
            self._initialized = True
    
    def observe(self, target_dict: Dict, callback: Union[Callable[[str, Any, Any], None], 
                                                        Callable[[str, Any, Any], asyncio.Future]]):
        """
        监听字典变化
        
        Args:
            target_dict: 要监听的目标字典
            callback: 回调函数，接收三个参数：操作类型('set', 'delete', 'update')、键、值
        """
        dict_id = id(target_dict)
        
        if dict_id not in self._watched_dicts:
            wrapped_dict = ObservableDict(target_dict, self._notify_callbacks)
            wrapped_id = id(wrapped_dict)
            
            self._watched_dicts[dict_id] = {
                'original': target_dict,
                'wrapped': wrapped_dict,
                'callbacks': set()
            }
            self._wrapped_to_original[wrapped_id] = dict_id
        
        self._watched_dicts[dict_id]['callbacks'].add(callback)
        return self._watched_dicts[dict_id]['wrapped']
    
    def disobserve(self, target_dict: Dict, callback: Union[Callable[[str, Any, Any], None], 
                                                           Callable[[str, Any, Any], asyncio.Future]] = None):
        """
        解除监听
        """
        dict_id = None
        if isinstance(target_dict, ObservableDict):
            wrapped_id = id(target_dict)
            if wrapped_id in self._wrapped_to_original:
                dict_id = self._wrapped_to_original[wrapped_id]
        else:
            dict_id = id(target_dict)
        
        if dict_id is not None and dict_id in self._watched_dicts:
            watched_info = self._watched_dicts[dict_id]
            
            if callback is None:
                watched_info['callbacks'].clear()
            else:
                watched_info['callbacks'].discard(callback)
            
            if len(watched_info['callbacks']) == 0:
                wrapped_id = id(watched_info['wrapped'])
                if wrapped_id in self._wrapped_to_original:
                    del self._wrapped_to_original[wrapped_id]
                del self._watched_dicts[dict_id]
    
    def _notify_callbacks(self, target_dict: Dict, operation: str, key: Any, value: Any):
        """通知所有回调函数"""
        dict_id = id(target_dict)
        
        if dict_id in self._watched_dicts:
            callbacks = self._watched_dicts[dict_id]['callbacks'].copy()
            
            for callback in callbacks:
                try:
                    if asyncio.iscoroutinefunction(callback):
                        # 异步回调处理
                        self._run_async_callback(callback, operation, key, value)
                    else:
                        # 同步回调：直接调用
                        callback(operation, key, value)
                        
                except Exception as e:
                    print(f"Callback error: {e}")
    
    def _run_async_callback(self, callback, operation, key, value):
        """运行异步回调函数"""
        try:
            # 尝试获取当前运行的事件循环
            loop = asyncio.get_running_loop()
            # 如果当前线程有事件循环，直接创建任务
            task = asyncio.create_task(callback(operation, key, value))
            
            # 添加异常处理
            def exception_handler(task):
                try:
                    task.result()
                except Exception as e:
                    print(f"Async callback error: {e}")
            
            task.add_done_callback(exception_handler)
            
        except RuntimeError:
            # 当前线程没有事件循环，需要在主线程的事件循环中运行
            # 创建一个新的事件循环并在后台运行
            def run_in_new_loop():
                new_loop = asyncio.new_event_loop()
                asyncio.set_event_loop(new_loop)
                try:
                    coro = callback(operation, key, value)
                    new_loop.run_until_complete(coro)
                except Exception as e:
                    print(f"Async callback error: {e}")
                finally:
                    new_loop.close()
            
            # 在新线程中运行事件循环
            thread = threading.Thread(target=run_in_new_loop, daemon=True)
            thread.start()


class ObservableDict(MutableMapping):
    """
    可观察的字典包装器
    """
    
    def __init__(self, original_dict: Dict, notify_func: Callable):
        self._data = original_dict
        self._notify = notify_func
        self._original_dict = original_dict
    
    def __getitem__(self, key):
        return self._data[key]
    
    def __setitem__(self, key, value):
        old_value = self._data.get(key, None) if key in self._data else None
        is_new = key not in self._data
        self._data[key] = value
        
        if is_new:
            self._notify(self._original_dict, 'set', key, value)
        else:
            self._notify(self._original_dict, 'update', key, (old_value, value))
    
    def __delitem__(self, key):
        value = self._data[key]
        del self._data[key]
        self._notify(self._original_dict, 'delete', key, value)
    
    def __iter__(self):
        return iter(self._data)
    
    def __len__(self):
        return len(self._data)
    
    def __repr__(self):
        return f"ObservableDict({self._data})"
    
    def update(self, *args, **kwargs):
        updates = {}
        
        if args:
            other = args[0]
            if hasattr(other, "items"):
                for key, value in other.items():
                    updates[key] = value
            else:
                for key, value in other:
                    updates[key] = value
        
        for key, value in kwargs.items():
            updates[key] = value
        
        for key, value in updates.items():
            self.__setitem__(key, value)

# 异步回调函数示例
async def async_callback(operation, key, value):
    print(f"异步回调 - 操作: {operation}, 键: {key}, 值: {value}")
    # 模拟异步操作
    await asyncio.sleep(0.01)
    print(f"异步回调完成 - {key}")

def sync_callback(operation, key, value):
    print(f"同步回调 - 操作: {operation}, 键: {key}, 值: {value}")

observer = DictObserver()

if __name__ == "__main__":
    import time
    import threading
    
    test_dict = {"name": "Alice", "age": 25}
    
    # 同时添加同步和异步回调
    observed_dict = observer.observe(test_dict, async_callback)
    observer.observe(test_dict, sync_callback)
    
    print("初始字典:", observed_dict)
    
    # 在主线程中测试
    observed_dict["city"] = "Beijing"
    observed_dict["age"] = 26
    del observed_dict["name"]
    
    # 在另一个线程中测试
    def thread_test():
        observed_dict["thread_key"] = "thread_value"
        observed_dict["age"] = 30
    
    thread = threading.Thread(target=thread_test)
    thread.start()
    thread.join()
    
    print("最终字典:", observed_dict)
    
    # 等待一段时间让异步回调完成
    time.sleep(0.5)
    
    print("测试完成")