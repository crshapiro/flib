# Linked list

The **list** class is a doubly-linked list that contain any class type. Intrinsic types must be returned using a wrapper function from the [standard library](stl.md).

## Constructor
*l* = **list_t**() generates an empty list.

## Methods
* **list_t%push_front**(*item*) Pushes *item* to the front of the list
* **list_t%push_back**(*item*) Pushes *item* to the back of the list
* **list_t%pop_front**(*item*) Pops *item* from the front of the list
* **list_t%pop_back**(*item*) Pops *item* from the back of the list
* **list_t%insert**(*item*, *n*) Inserts *item* into the *n*th location of the list
* **list_t%remove**(*n*) Removes and returns the *n*th item from the list
* **list_t%size**() Returns the number of items in the list