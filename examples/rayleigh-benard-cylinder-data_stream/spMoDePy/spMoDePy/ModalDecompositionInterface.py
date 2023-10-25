import abc

class ModalDecompositionInterface(metaclass=abc.ABCMeta):
    @classmethod
    def __subclasshook__(cls, subclass):
        return (hasattr(subclass, 'firststep') and
                callable(subclass.firststep) and
                hasattr(subclass, 'loadbuffer1') and
                callable(subclass.loadbuffer1) and
                hasattr(subclass, 'loadbuffer2') and
                callable(subclass.loadbuffer2) and
                hasattr(subclass, 'update') and
                callable(subclass.update) and
                hasattr(subclass, 'step') and
                callable(subclass.step) and
                hasattr(subclass, 'poststream') and
                callable(subclass.poststream) or
                NotImplemented)

    @abc.abstractmethod
    def firststep(self, params,case,iterator, io_controller,comm):
        """Describe later"""
        raise NotImplementedError

    @abc.abstractmethod
    def load_buffer1(self, params, case, iterator, io_controller, comm):
        """Describe later"""
        raise NotImplementedError
    
    @abc.abstractmethod
    def update_from_buffer1(self, params, case, iterator, io_controller, comm):
        """Describe later"""
        raise NotImplementedError
    
    @abc.abstractmethod
    def load_buffer2(self, *args, **kwargs):
        """Describe later"""
        raise NotImplementedError 
    
    @abc.abstractmethod
    def update_from_buffer2(self, *args, **kwargs):
        """Describe later"""
        raise NotImplementedError
    
    @abc.abstractmethod
    def poststream(self, params, case, iterator, io_controller, comm):
        """Describe later"""
        raise NotImplementedError
