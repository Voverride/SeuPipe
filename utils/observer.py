from typing import Dict, Callable, Any, Union, Optional
from collections.abc import MutableMapping
import asyncio
import threading


class DictObserver:
    """
    字典变化观察器 - 单例模式（支持异步回调和额外参数传递）
    """
    _instance = None
    _initialized = False
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(DictObserver, cls).__new__(cls)
        return cls._instance
    
    def __init__(self):
        if not self._initialized:
            # 使用字典ID作为键，避免弱引用问题
            self._watched_dicts = {}  # {id(original_dict): info}
            self._wrapped_to_original_id = {}  # {id(wrapped_dict): original_dict_id}
            self._initialized = True
    
    def observe(self, 
                target_dict: Dict, 
                callback: Union[Callable, asyncio.Future],
                **kwargs) -> MutableMapping:
        """
        监听字典变化，支持传递额外参数
        
        Args:
            target_dict: 要监听的目标字典
            callback: 回调函数，接收参数: (operation, key, value, **kwargs)
            **kwargs: 任意额外参数，将在回调时传递给callback
            
        Returns:
            包装后的可观察字典
        """
        # 如果已经是被包装的字典，获取原始字典
        if isinstance(target_dict, ObservableDict):
            target_dict = target_dict._original_dict
        
        original_dict_id = id(target_dict)
        
        if original_dict_id not in self._watched_dicts:
            wrapped_dict = ObservableDict(target_dict, self._notify_callbacks)
            wrapped_id = id(wrapped_dict)
            
            # 存储包装字典和原始字典ID的映射
            self._wrapped_to_original_id[wrapped_id] = original_dict_id
            
            self._watched_dicts[original_dict_id] = {
                'wrapped': wrapped_dict,
                'callbacks': []  # 存储元组: (callback, kwargs)
            }
        else:
            wrapped_dict = self._watched_dicts[original_dict_id]['wrapped']
        
        # 检查是否已存在相同回调（避免重复注册）
        callbacks = self._watched_dicts[original_dict_id]['callbacks']
        callback_exists = False
        
        for cb, cb_kwargs in callbacks:
            if cb is callback and cb_kwargs == kwargs:
                callback_exists = True
                break
        
        if not callback_exists:
            callbacks.append((callback, kwargs))
        
        return wrapped_dict
    
    def disobserve(self, 
                  target_dict: Dict, 
                  callback: Optional[Callable] = None, 
                  **kwargs) -> None:
        """
        解除监听
        
        Args:
            target_dict: 要解除监听的字典（原始字典或包装字典）
            callback: 要移除的回调函数。如果为None，则移除所有回调
            **kwargs: 与注册时相同的额外参数。如果提供，则只移除匹配参数的回调
        """
        # 确定原始字典ID
        if isinstance(target_dict, ObservableDict):
            wrapped_id = id(target_dict)
            if wrapped_id in self._wrapped_to_original_id:
                original_dict_id = self._wrapped_to_original_id[wrapped_id]
            else:
                return
        else:
            original_dict_id = id(target_dict)
        
        if original_dict_id not in self._watched_dicts:
            return
        
        watched_info = self._watched_dicts[original_dict_id]
        
        if callback is None:
            # 移除所有回调
            watched_info['callbacks'].clear()
        else:
            # 筛选保留不匹配的回调
            new_callbacks = []
            for cb, cb_kwargs in watched_info['callbacks']:
                # 保留不匹配的回调
                if not (cb is callback and (not kwargs or cb_kwargs == kwargs)):
                    new_callbacks.append((cb, cb_kwargs))
            watched_info['callbacks'] = new_callbacks
        
        # 如果没有回调了，清理包装
        if not watched_info['callbacks']:
            wrapped_dict = watched_info['wrapped']
            wrapped_id = id(wrapped_dict)
            if wrapped_id in self._wrapped_to_original_id:
                del self._wrapped_to_original_id[wrapped_id]
            del self._watched_dicts[original_dict_id]
    
    def _notify_callbacks(self, 
                         target_dict: Dict, 
                         operation: str, 
                         key: Any, 
                         value: Any) -> None:
        """通知所有回调函数"""
        original_dict_id = id(target_dict)
        
        if original_dict_id not in self._watched_dicts:
            return
        
        callbacks = self._watched_dicts[original_dict_id]['callbacks'].copy()
        
        for callback, kwargs in callbacks:
            try:
                # 检查是否是异步函数
                is_async = asyncio.iscoroutinefunction(callback) or (
                    hasattr(callback, '__call__') and 
                    asyncio.iscoroutinefunction(callback.__call__)
                )
                
                if is_async:
                    # 异步回调处理
                    self._run_async_callback(callback, operation, key, value, **kwargs)
                else:
                    # 同步回调
                    callback(operation, key, value, **kwargs)
                    
            except Exception as e:
                print(f"Callback error: {e}")
    
    def _run_async_callback(self, callback, operation, key, value, **kwargs):
        """运行异步回调函数"""
        # 创建异步任务
        async def run_callback():
            try:
                await callback(operation, key, value, **kwargs)
            except Exception as e:
                print(f"Async callback error: {e}")
        
        try:
            # 尝试获取当前事件循环
            loop = asyncio.get_running_loop()
            # 创建任务
            task = loop.create_task(run_callback())
            # 添加异常处理
            task.add_done_callback(lambda t: t.exception() and 
                                  print(f"Async callback task failed: {t.exception()}"))
        except RuntimeError:
            # 没有运行中的事件循环，创建新线程运行
            def run_in_new_loop():
                new_loop = asyncio.new_event_loop()
                asyncio.set_event_loop(new_loop)
                try:
                    new_loop.run_until_complete(run_callback())
                except Exception as e:
                    print(f"Async callback error in new loop: {e}")
                finally:
                    new_loop.close()
            
            thread = threading.Thread(target=run_in_new_loop, daemon=True)
            thread.start()


class ObservableDict(MutableMapping):
    """
    可观察的字典包装器
    """
    
    __slots__ = ('_data', '_notify', '_original_dict')
    
    def __init__(self, original_dict: Dict, notify_func: Callable):
        self._data = original_dict
        self._notify = notify_func
        self._original_dict = original_dict
    
    def __getitem__(self, key):
        return self._data[key]
    
    def __setitem__(self, key, value):
        old_value = self._data.get(key) if key in self._data else None
        is_new = key not in self._data
        self._data[key] = value
        
        if is_new:
            self._notify(self._original_dict, 'set', key, value)
        else:
            self._notify(self._original_dict, 'update', key, (old_value, value))
    
    def __delitem__(self, key):
        if key not in self._data:
            raise KeyError(key)
        value = self._data.pop(key)
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
                updates.update(other.items())
            elif hasattr(other, "__iter__"):
                for item in other:
                    if len(item) == 2:
                        k, v = item
                        updates[k] = v
                    else:
                        raise ValueError("update() takes either a dict-like object or an iterable of (key, value) pairs")
        
        updates.update(kwargs)
        
        for key, value in updates.items():
            self[key] = value

observer = DictObserver()



# 使用示例
# 异步回调函数示例
async def async_callback(operation, key, value, taskName=None, slice=None, clipName=None):
    print(f"Async callback - Task: {taskName}, Slice: {slice}, Clip: {clipName}")
    print(f"  Operation: {operation}, Key: {key}, Value: {value}")
    await asyncio.sleep(0.1)
    print(f"  Async callback completed for {key}")

# 同步回调函数示例
def sync_callback(operation, key, value, extra_info=None):
    print(f"Sync callback - Extra: {extra_info}")
    print(f"  Operation: {operation}, Key: {key}, Value: {value}")

if __name__ == "__main__":
    import asyncio
    import time
    import threading
    status = {"initial": "value"}
    
    # 监听字典并传递额外参数
    observed_status = observer.observe(
        status,
        async_callback,
        taskName="render_task",
        slice="main_slice",
        clipName="intro_clip"
    )
    
    observer.observe(
        status,
        sync_callback,
        extra_info="important_data"
    )
    
    print("初始状态:", dict(observed_status))
    
    # 测试修改
    observed_status["new_key"] = "new_value"  # 触发'set'
    observed_status["initial"] = "updated"    # 触发'update'
    
    # 在子线程中修改
    def thread_test():
        observed_status["thread_key"] = "from_thread"
        del observed_status["initial"]  # 触发'delete'
    
    thread = threading.Thread(target=thread_test)
    thread.start()
    thread.join()
    
    # 等待异步回调完成
    time.sleep(0.5)
    
    # 解除特定回调
    observer.disobserve(
        observed_status,
        async_callback,
        taskName="render_task",
        slice="main_slice",
        clipName="intro_clip"
    )
    
    # 验证解除后不再触发
    observed_status["after_removal"] = "should_not_trigger_async"
    
    print("最终状态:", dict(observed_status))
    print("测试完成")